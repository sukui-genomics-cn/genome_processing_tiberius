import argparse
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
    def __init__(self, index_file, shuffle=False, max_maps: int = 8, mode="r", seq_key="seq", anno_key="anno"):
        """
        Args:
            index_file: Preprocessed chunk index file
            shuffle: Whether to shuffle the data order
        """
        self.chunks = []
        self.file_handles = {}  # Mapping from file path to file descriptor
        self.lock = torch.multiprocessing.Lock()  # Used for multi-process synchronization
        self.max_maps = max_maps
        self.mode = mode
        self.seq_key = seq_key
        self.anno_key = anno_key
        self.file_cache = OrderedDict() 
        
        # Load index
        self.chunks = self._load_chunk_file(index_file)
        
        if shuffle:
            random.shuffle(self.chunks)
            self.chunks = sorted(self.chunks, key=lambda x: (x[0]))  # Sort by chromosome name
        
        # Register cleanup function
        atexit.register(self.cleanup)
    
    def _load_chunk_file(self, chunk_file):
        """Load chunk.txt file and return structured chunk information"""
        chunks = []
        with open(chunk_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    chunks.append([parts[0], int(parts[1]), int(parts[2])])
        logger.info(f"Loaded {len(chunks)} chunks from {chunk_file}")    
        return chunks

    def _get_mmap(self, file_path, anno_path=None):
        """Get or create memory mapping (thread-safe)"""
        with self.lock:
            if file_path not in self.file_cache:
                if not os.path.exists(file_path):
                    raise FileNotFoundError(f"SEQ file not found: {file_path}")
                
                # Open file in read mode and build mmap
                self.file_cache[file_path] = h5py.File(file_path, self.mode, driver='core' if self.mode == 'r' else None)
                self.file_cache.move_to_end(file_path)  # Move to end to maintain order

            if len(self.file_cache) > self.max_maps:
                print(f"Max maps reached, removing oldest mmap: {self.max_maps}")

            # Check mmap size
            if len(self.file_cache) > self.max_maps:
                self._remove_oldest_mmap()
                
            return self.file_cache[file_path]
        
    def _remove_oldest_mmap(self):
        """Close the oldest mmap and remove it from the dictionary"""
        if not self.file_cache:
            return

        oldest_path, oldest_mmap = next(iter(self.file_cache.items()))
        if oldest_path in self.file_cache:
            # Close mapping and file descriptor
            oldest_mmap.close()
            del self.file_cache[oldest_path]

        print(f"Released oldest memory mapping: {oldest_path}")
        print(f"Current retained memory mappings: {self.file_cache.keys()}")
    
    def _encode_dna(self, sequence):
        """Convert DNA sequence to one-hot encoding"""
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
        """Clean up resources"""
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
    """Get DataLoader for DNA sequences"""
    # Set fork method for shared memory
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
        # multiprocessing_context=torch.multiprocessing.get_context('spawn')  # Use fork instead of spawn
    )
    
    return dataloader

if __name__ == "__main__":
    # Example usage
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset_file', help='Path to the dataset file')
    
    args = parser.parse_args()
    dataloader = get_dna_dataloader(args.dataset_file, batch_size=32, num_workers=0, shuffle=True)
    
    for batch in dataloader:
        print(batch["input_seq"].shape, batch["anno"].shape if batch["anno"] is not None else "No annotations")
        break  # Test the first batch only