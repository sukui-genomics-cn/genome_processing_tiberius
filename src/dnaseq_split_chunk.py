import argparse
import glob
import os
import random
from Bio import SeqIO
import logging

import h5py
from tqdm import tqdm

logger = logging.getLogger(__name__)

class SequencePreprocessor:
    """预处理DNA序列文件并生成索引的类"""
    
    CHROMOSOME_KEYWORDS = ['chr', 'chromosome', 'NC_', 'NT_']
    EXCLUDE_KEYWORDS = ['scaffold', 'contig', 'random', 'Un']
    SUPPORTED_FASTA_EXTENSIONS = (".h5")
    
    def __init__(self, input_dir, output_dir, chunk_size=500000):
        """
        初始化预处理类
        
        Args:
            input_dir: 输入目录路径
            output_dir: 输出目录路径
            chunk_size: 分块大小，默认为50000
        """
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.chunk_size = chunk_size
        
        # 确保输入目录存在
        self.h5_input_dir = os.path.join(input_dir, "chunk_chr.all")
        self.chr_files = self.find_chr_files()
        # 输出文件路径
        self.index_file = os.path.join(output_dir, "chunk_index.txt")
        self.meta_file = os.path.join(output_dir, "metadata.txt")
    
    def find_chr_files(self, format="*.h5"):
        """查找所有染色体文件"""
        chr_files = glob.glob(os.path.join(self.h5_input_dir, format))
        logger.info(f"Found {len(chr_files)} chromosome files")
        return chr_files

    @classmethod
    def is_chr(cls, record_id):
        """
        判断给定的记录ID是否属于染色体
        
        Args:
            record_id: 记录ID字符串
            
        Returns:
            bool: 如果是染色体返回True，否则返回False
        """
        # 规则1：包含染色体关键词
        is_chr = any(kw in record_id for kw in cls.CHROMOSOME_KEYWORDS)
        # 规则2：不包含排除关键词
        is_clean = not any(ex in record_id for ex in cls.EXCLUDE_KEYWORDS)
        return is_chr and is_clean
    
    def _process_h5_file(self, h5_file):
        """处理文本文件(.bin或.txt)"""
        if not os.path.exists(h5_file):
            logger.warning(f"h5 file does not exist: {h5_file}")
            return
        
        logger.info(f"Processing text file: {h5_file}")
        chunk_list = []
        meta_list = []
        with h5py.File(h5_file, 'r') as input_f:
            seq_len = input_f['seq'].shape[0]
            if seq_len < self.min_seq_len:
                logger.info(f"Sequence length {seq_len} of {h5_file} is too short, skipping.")
            else:
                num_chunks = seq_len // self.chunk_size
                
                # 记录索引
                for i in range(num_chunks):
                    start = i * self.chunk_size
                    end = start + self.chunk_size
                    chunk_list.append([h5_file, start, end])
                
                meta_list = [h5_file, seq_len, self.chunk_size, len(chunk_list)]
        return chunk_list, meta_list

    
    def process(self, chr_shuffle=True, species_shuffle=True):
        """执行预处理流程"""
        os.makedirs(self.output_dir, exist_ok=True)
        chunks = []
        meatas = []
            
        for chr_file in tqdm(self.chr_files, desc="Processing files", total=len(self.chr_files)):
            chunk_list, meta_list = self._process_h5_file(chr_file)
            if chunk_list:
                if chr_shuffle:
                    random.shuffle(chunk_list)
                
                chunks.append(chunk_list)
                if meta_list:
                    meatas.extend(meta_list)
        if species_shuffle:
            random.shuffle(chunks)

        # 写入索引文件
        with open(self.index_file, 'w') as idx_f:
            for chunk in chunks:
                for item in chunk:
                    idx_f.write(f"{item[0]}\t{item[1]}\t{item[2]}\n")
        # 写入元数据文件
        with open(self.meta_file, 'w') as meta_f:
            for meta in meatas:
                meta_f.write(f"{meta[0]}\t{meta[1]}\t{meta[2]}\t{meta[3]}\n")
        
        logger.info(f"Index file created at: {self.index_file} and {self.meta_file}")
                

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", type=str, required=True)
    parser.add_argument("--output_dir", type=str, required=True)
    parser.add_argument("--chunk_size", type=int, default=50000)
    args = parser.parse_args()

    preprocessor = SequencePreprocessor(
        input_dir=args.input_dir, 
        output_dir=args.output_dir,
        chunk_size=args.chunk_size
        )
    preprocessor.process()


if __name__ == "__main__":
    # 使用示例
    main()