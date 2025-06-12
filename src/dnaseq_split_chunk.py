import argparse
import copy
import glob
import os
import random
from Bio import SeqIO
import logging

import h5py
from tqdm import tqdm

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class SequencePreprocessor:
    """预处理DNA序列文件并生成索引的类"""
    
    CHROMOSOME_KEYWORDS = ['chr', 'chromosome', 'NC_', 'NT_']
    EXCLUDE_KEYWORDS = ['scaffold', 'contig', 'random', 'Un']
    SUPPORTED_FASTA_EXTENSIONS = (".h5")
    
    def __init__(self, input_dir:str, output_dir:str, sub_dir="**", recursive:bool=True, 
                 filter_file_path:str=None, name:str="*.h5", chunk_size:int=500000, 
                 min_seq_len:int=1000000, seed=42):
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
        self.min_seq_len = min_seq_len
        
        # 确保输入目录存在

        self.chr_files = self.find_chr_files(
            dest_path=self.input_dir,  # 输入目录
            sub_dir=sub_dir,
            recursive=recursive,  # 是否递归查找子目录
            filter_file_path=filter_file_path,  # 可以指定过滤文件路径
            name=name  # 查找所有的h5文件
        )
        # 输出文件路径
        self.chunks = []
        self.meatas = []

        random.seed(seed)
        
    
    def find_chr_files(self, dest_path, sub_dir,recursive , filter_file_path, name):
        """查找所有染色体文件"""
        # chr_files = glob.glob(os.path.join(self.h5_input_dir, format), recursive=True)
        # logger.info(f"Found {len(chr_files)} chromosome files")


        if filter_file_path is not None and os.path.exists(filter_file_path):
            with open(filter_file_path, "r", encoding="utf8") as f:
                filter_files = f.readlines()
            filter_path_root = os.path.dirname(filter_file_path)
            filter_files = [os.path.join(filter_path_root, file.strip()) for file in filter_files]
        else:
            filter_files = []

        assert os.path.exists(dest_path), f"{dest_path} is not exits"
        file_dirs = []

        if recursive:
            for sub_file_name in os.listdir(dest_path):
                sub_file_dir = os.path.join(dest_path, sub_file_name)
                if os.path.isdir(sub_file_dir) and sub_file_dir not in filter_files:
                    file_dirs += glob.glob(os.path.join(dest_path, sub_file_dir, sub_dir, name), recursive=recursive)
                else:
                    logger.info(f"skip file: {sub_file_dir}")
        else:
            if sub_dir is not None:
                logger.info(f"sub_dir: {sub_dir}, not support when recursive is False")
            file_dirs += glob.glob(os.path.join(dest_path, name), recursive=recursive)
            logger.info(f"find files: {len(file_dirs)}")

        return file_dirs

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

    
    def process(self):
        """执行预处理流程"""
        os.makedirs(self.output_dir, exist_ok=True)
        
        for chr_file in tqdm(self.chr_files, desc="Processing files", total=len(self.chr_files)):
            chunk_list, meta_list = self._process_h5_file(chr_file)
            if chunk_list:
                self.chunks.append(chunk_list)
            if meta_list:
                self.meatas.append(meta_list)
        self.write_index(self.chunks, self.meatas, self.output_dir, name_postfix="all")
        
    
    def write_index(self, chunks:list,meatas:list, output_dir:str, name_postfix:str="all"):
        """将处理后的数据写入索引文件和元数据文件"""
        # 写入索引文件
        index_file = os.path.join(output_dir, f"chunks_{name_postfix}.txt")
        with open(index_file, 'w') as idx_f:
            for chunk in chunks:
                for item in chunk:
                    idx_f.write(f"{item[0]}\t{item[1]}\t{item[2]}\n")
        # 写入元数据文件
        meta_file = os.path.join(output_dir, f"meta_{name_postfix}.txt")
        with open(meta_file, 'w') as meta_f:
            for meta in meatas:
                meta_f.write(f"{meta[0]}\t{meta[1]}\t{meta[2]}\t{meta[3]}\n")
        
        logger.info(f"Index file created at: {index_file} and {meta_file}")

    def split_train_val_test(self, dataset_name:str, train_ratio=0.8, val_ratio=0.1, test_ratio=0.1, shuffle=True):
        """将数据集划分为训练集、验证集和测试集"""
        output_dir = os.path.join(self.output_dir, dataset_name)
        os.makedirs(output_dir, exist_ok=True)

        # split chunks into train, val, test
        chunks = copy.deepcopy(self.chunks)
        assert len(chunks) == len(self.meatas), "Chunks and metadatas must have the same length"
        indexs = list(range(len(chunks)))
        if shuffle:
            random.shuffle(indexs)
        
        chunks_train, metas_train = [], []
        chunks_val, metas_val = [], []
        chunks_test, metas_test = [], []
        for i, (indexs) in enumerate(indexs):
              # set random seed for reproducibility
            rate = random.random()
            if rate < train_ratio:
                chunk = chunks[indexs]
                meta = self.meatas[indexs]
                if shuffle:
                    random.shuffle(chunk)
                chunks_train.append(chunk)
                metas_train.append(meta)
            elif rate > train_ratio + val_ratio:
                chunk = chunks[indexs]
                meta = self.meatas[indexs]
                if shuffle:
                    random.shuffle(chunk)
                chunks_test.append(chunk)
                metas_test.append(meta)
            else:
                chunk = chunks[indexs]
                meta = self.meatas[indexs]
                if shuffle:
                    random.shuffle(chunk)
                chunks_val.append(chunk)
                metas_val.append(meta)

        # write train, val, test index files
        self.write_index(chunks_train, metas_train, output_dir, name_postfix="train")
        self.write_index(chunks_val, metas_val, output_dir, name_postfix="val")
        self.write_index(chunks_test, metas_test, output_dir, name_postfix="test")
        logger.info(f"Train, val, test index files created at: {output_dir}")
        logger.info(f"Train: {len(chunks_train)}, Val: {len(chunks_val)}, Test: {len(chunks_test)}")                               
            

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", type=str, required=True)
    parser.add_argument("--sub_dir", type=str, default="**", help="Subdirectory to search for files")
    parser.add_argument("--recursive", action='store_true', help="Whether to search recursively in subdirectories")
    parser.add_argument("--filter_file_path", type=str, default=None, help="Path to a file containing paths to filter out")
    parser.add_argument("--name", type=str, default="*.h5", help="File name pattern to match")
    parser.add_argument("--output_dir", type=str, required=True)
    parser.add_argument("--chunk_size", type=int, default=50000)
    parser.add_argument("--min_seq_len", type=int, default=1000000, help="Minimum sequence length to process")
    parser.add_argument("--train_ratio", type=float, default=0.8, help="Ratio of training set")
    parser.add_argument("--val_ratio", type=float, default=0.1, help="Ratio of validation set")
    parser.add_argument("--test_ratio", type=float, default=0.1, help="Ratio of test set")
    parser.add_argument("--shuffle", action='store_true', help="Shuffle the dataset before splitting")
    parser.add_argument("--dataset_name", type=str, default="dnaseq_split", help="Name of the dataset for splitting")

    args = parser.parse_args()
    logger.setLevel(logging.INFO)
    logger.info(f"args: {args}")
    # 创建输出目录

    preprocessor = SequencePreprocessor(
        input_dir=args.input_dir,
        sub_dir=args.sub_dir,
        recursive=args.recursive,
        filter_file_path=args.filter_file_path,
        name=args.name,
        output_dir=args.output_dir,
        chunk_size=args.chunk_size
        )
    preprocessor.process()
    preprocessor.split_train_val_test(
        dataset_name=args.dataset_name, 
        train_ratio=args.train_ratio, 
        val_ratio=args.val_ratio, 
        test_ratio=args.test_ratio, 
        shuffle=args.shuffle
        )


if __name__ == "__main__":
    # 使用示例
    main()