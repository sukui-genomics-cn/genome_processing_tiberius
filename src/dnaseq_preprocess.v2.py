import argparse
import os
from Bio import SeqIO
import logging

from tqdm import tqdm

logger = logging.getLogger(__name__)

class SequencePreprocessor:
    """预处理DNA序列文件并生成索引的类"""
    
    CHROMOSOME_KEYWORDS = ['chr', 'chromosome', 'NC_', 'NT_']
    EXCLUDE_KEYWORDS = ['scaffold', 'contig', 'random', 'Un']
    SUPPORTED_FASTA_EXTENSIONS = ('.txt', '.bin', '.fasta', '.fna')
    
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
        self.fasta_input_dir = os.path.join(input_dir, "chunk_chr.all")
        self.anno_input_dir = os.path.join(input_dir, "chunk_chr.all")
        self._validate_input_dirs()
        
        # 输出文件路径
        self.index_file = os.path.join(output_dir, "chunk_index.txt")
        self.meta_file = os.path.join(output_dir, "metadata.txt")
        
    def _validate_input_dirs(self):
        """验证输入目录是否存在"""
        if not os.path.exists(self.fasta_input_dir):
            raise FileNotFoundError(f"FASTA input directory does not exist: {self.fasta_input_dir}")
        if not os.path.exists(self.anno_input_dir):
            raise FileNotFoundError(f"anno_tiberius input directory does not exist: {self.anno_input_dir}")
    
    @classmethod
    def is_chromosome(cls, record_id):
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
    
    def _process_text_file(self, filepath, anno_path, idx_f, meta_f):
        """处理文本文件(.bin或.txt)"""
        if not os.path.exists(anno_path):
            logger.warning(f"PKL file does not exist: {anno_path}")
            return
        
        logger.info(f"Processing text file: {filepath}")
        with open(filepath, 'r') as input_f:
            seq = input_f.read().strip()
            seq_len = len(seq)
            num_chunks = seq_len // self.chunk_size
            
            # 记录索引
            for i in range(num_chunks):
                start = i * self.chunk_size
                end = start + self.chunk_size
                idx_f.write(f"{filepath}\t{anno_path}\t{start}\t{end}\n")
            
            meta_f.write(f"{filepath}\t{anno_path}\t{seq_len}\n")
    
    def _process_fasta_file(self, filepath, idx_f, meta_f):
        """处理FASTA文件(.fasta或.fna)"""
        logger.info(f"Processing FASTA file: {filepath}")
        filename = os.path.basename(filepath)
        output_file = os.path.join(self.output_dir, f"{os.path.splitext(filename)[0]}.bin")
        
        with open(output_file, 'wb') as out_f:
            for record in SeqIO.parse(filepath, "fasta"):
                if not self.is_chromosome(record.id):
                    logger.info(f"Skipping non-chromosome record: {record.id} in {filename}")
                    continue
                    
                seq = str(record.seq)
                seq_len = len(seq)
                num_chunks = seq_len // self.chunk_size
                
                # 写入二进制数据
                out_f.write(seq.encode('ascii'))
                
                # 记录索引
                for i in range(num_chunks):
                    start = i * self.chunk_size
                    end = start + self.chunk_size
                    idx_f.write(f"{output_file}\t{output_file}\t{start}\t{end}\n")
                
                meta_f.write(f"{filename}\t{filename}\t{seq_len}\n")
    
    def process(self):
        """执行预处理流程"""
        os.makedirs(self.output_dir, exist_ok=True)
        
        with open(self.index_file, 'w') as idx_f, open(self.meta_file, 'w') as meta_f:
            meta_f.write(f"chunk_size\t{self.chunk_size}\n")
            
            for filename in tqdm(os.listdir(self.fasta_input_dir), desc="Processing files", total=len(os.listdir(self.fasta_input_dir))):
                print(f"Processing file: {filename}")
                if not filename.endswith(self.SUPPORTED_FASTA_EXTENSIONS):
                    logger.info(f"Skipping unsupported file type: {filename}")
                    continue
                    
                filepath = os.path.join(self.fasta_input_dir, filename)
                
                if filename.endswith(('.bin', '.txt')):
                    anno_name = filename.replace(".bin", ".npy")
                    anno_path = os.path.join(self.anno_input_dir, anno_name)
                    self._process_text_file(filepath, anno_path, idx_f, meta_f)
                elif filename.endswith(('.fna', '.fasta')):
                    self._process_fasta_file(filepath, idx_f, meta_f)

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