import os
import mmap
import random
import torch
import numpy as np
from torch.utils.data import Dataset, DataLoader
from torch.utils.data.distributed import DistributedSampler
import atexit
from collections import OrderedDict
import logging
import h5py

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class DNAH5Dataset(Dataset):
    def __init__(self, index_file, shuffle=True, max_maps:int=8, mode="r", seq_key="seq", anno_key="anno"):
        """
        Args:
            index_file: 预处理生成的chunk索引文件
            shuffle: 是否打乱数据顺序
        """
        self.chunks = []
        self.file_handles = {}  # 文件路径到文件描述符的映射
        self.lock = torch.multiprocessing.Lock()  # 用于多进程同步
        self.max_maps = max_maps
        self.mode = mode
        self.seq_key = seq_key
        self.anno_key = anno_key
        self.file_cache = OrderedDict() 
        
        # 加载索引
        self.chunks = self._load_chunk_file(index_file)
        
        if shuffle:
            random.shuffle(self.chunks)
        
        # 注册清理函数
        atexit.register(self.cleanup)
    
    def _load_chunk_file(self, chunk_file):
        """加载chunk.txt文件，返回结构化的chunk信息"""
        chunks = []
        with open(chunk_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    chunks.append([parts[0], int(parts[1]), int(parts[2])])
        logger.info(f"Loaded {len(chunks)} chunks from {chunk_file}")    
        return chunks

    def _get_mmap(self, file_path, anno_path=None):
        """获取或创建内存映射（线程安全）"""
        with self.lock:
            if file_path not in self.file_cache:
                if not os.path.exists(file_path):
                    raise FileNotFoundError(f"SEQ file is not found: {file_path}")
                
                # open file by r and build mmap
                self.file_cache[file_path] = h5py.File(file_path, self.mode, driver='core' if self.mode == 'r' else None)
                self.file_cache.move_to_end(file_path) # mv to end to keep order

            if len(self.file_cache) > self.max_maps:
                print(f"max maps reached, will remove oldest mmap: {self.max_maps}")

            # check mmap size
            if len(self.file_cache) > self.max_maps:
                self._remove_oldest_mmap()
                
            return self.file_cache[file_path]
        
    def _remove_oldest_mmap(self):
        """close the oldest mmap and remove it from the dictionary"""
        if not self.file_cache:
            return

        oldest_path, oldest_mmap = next(iter(self.file_cache.items()))
        if oldest_path in self.file_cache:
            # 关闭映射和文件描述符
            oldest_mmap.close()
            del self.file_cache[oldest_path]

        print(f"已释放最旧的内存映射: {oldest_path}")
        print(f"当前保留的内存映射: {self.file_cache.keys()}")
    
    def _encode_dna(self, sequence):
        """将DNA序列转换为one-hot编码"""
        mapping = {
            'A': [1, 0, 0, 0],
            'T': [0, 1, 0, 0],
            'C': [0, 0, 1, 0],
            'G': [0, 0, 0, 1],
            'N': [0.25, 0.25, 0.25, 0.25]
        }
        
        encoded = np.zeros((len(sequence), 4), dtype=np.float32)
        for i, base in enumerate(sequence):
            encoded[i] = mapping.get(base, mapping['N'])
        
        return torch.from_numpy(encoded)
    
    def __len__(self):
        return len(self.chunks)
    
    def __getitem__(self, idx):
        file_path, start, end = self.chunks[idx]
        
        try:
            data = self._get_mmap(file_path)
            seq = data[self.seq_key][start:end].tobytes().decode('ascii')
            anno = data[self.anno_key][start:end] if self.anno_key in data else None
            seq = self._encode_dna(seq)
            
            return {"input_seq": seq, "anno": anno}
        
        except Exception as e:
            print(f"Error processing {file_path}[{start}:{end}]: {str(e)}")
            raise
    
    def cleanup(self):
        """清理资源"""
        with self.lock:
            for file_path, mm in self.file_cache.items():
                try:
                    mm.close()
                except:
                    pass
            
            self.file_cache.clear()
    
    def __del__(self):
        self.cleanup()

def get_dna_dataloader(index_file, batch_size=32, num_workers=4, shuffle=True, distributed=False):
    """获取DNA序列的DataLoader"""
    # 设置共享内存的fork方式
    torch.multiprocessing.set_sharing_strategy('file_system')
    
    dataset = DNAH5Dataset(index_file, shuffle=shuffle)
    
    sampler = None
    if distributed:
        sampler = DistributedSampler(dataset, shuffle=shuffle)
    
    dataloader = DataLoader(
        dataset,
        batch_size=batch_size,
        sampler=sampler,
        num_workers=num_workers,
        pin_memory=True,
        persistent_workers=num_workers > 0,
        drop_last=True,
        # multiprocessing_context=torch.multiprocessing.get_context('spawn')  # 使用fork而非spawn
    )
    
    return dataloader

if __name__ == "__main__":
    # Example usage
    index_file = "/home/nvme01/sukui/01.data/T2T/t2t_chr_50004/chunk_index.txt"
    dataloader = get_dna_dataloader(index_file, batch_size=32, num_workers=0, shuffle=True)
    
    for batch in dataloader:
        print(batch["input_seq"].shape, batch["anno"].shape if batch["anno"] is not None else "No annotations")
        break  # Just to test the first batch