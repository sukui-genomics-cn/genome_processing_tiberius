import h5py
import os
import numpy as np

def test_h5():
    save_dir="/home/nvme01/data/T2T/chunk_chr.h5"
    h5_file = "/home/nvme01/data/T2T/chunk_chr.h5/T2T_NC_060925.1_backward.h5"


    # with h5py.File(os.path.join(save_dir, 'test.h5'), 'w') as f:
    #     dset = f.create_dataset('test', (10,), dtype='i')
    #     dset[:] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    
    with h5py.File(h5_file, 'r') as f:
        print(f"key: {f.keys()}")
        print(f"anno shape: {f['anno'].shape}")
        print(f"anno: {f['anno'][:10]}")
        print(f"seq shape: {f['seq'].shape}")
        print(f"seq: {f['seq'][:10]}")
        print(f"seq: {f['seq'][:10].tobytes().decode('ascii')}")


if __name__ == '__main__':
    test_h5()