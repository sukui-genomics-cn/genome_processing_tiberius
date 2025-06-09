import argparse
import glob
import os
import random
import logging
import numpy as np

import tqdm

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)



def parse_args():
    parse = argparse.ArgumentParser(description='Probiotic data analysis')
    parse.add_argument("-dp", '--dest_path', required=True, type=str, help='negative_path')
    parse.add_argument("-sp", '--save_path', required=True, type=str, help='save path')
    parse.add_argument("-tr", '--train_ratio', required=False, default=0.95, type=float, help='train data ratio')
    parse.add_argument("-vr", '--val_ratio', required=False, default=0.01, type=float, help='val data ratio')
    parse.add_argument("-ttr", '--test_ratio', required=False, default=0.04, type=float, help='test data ratio')
    parse.add_argument("-fa", '--file_name', required=False, default="*.pkl", type=str, help='file name')
    parse.add_argument("-rc", '--recursive', required=False, default=False, type=bool, help='recursive')
    parse.add_argument("-ffp", '--filter_file_path', required=False, default=None, type=str,
                       help='file of filter file name')
    args = parse.parse_args()
    logger.info(f"args: {args}")
    return args


class T2TDataProcess:
    def __init__(self, seed=42):
        self.seed = seed
        self._set_random_seed(self.seed)

    def split_pretrain_dta(
            self,
            dest_path: str = "/home/share/huadjyin/home/sunhaotong/trash/Growth-promoting_Bacteria/",
            save_path: str = "/home/share/huadjyin/home/sunhaotong/trash",
            train_ratio: float = 0.9,
            val_ratio: float = 0.05,
            test_ratio: float = 0.05,
            name: str = "*.pkl",
            recursive: bool = True,
            filter_file_path: str = None
    ):
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
                    file_dirs += glob.glob(os.path.join(dest_path, sub_file_dir, "**", name), recursive=recursive)
                else:
                    logger.info(f"skip file: {sub_file_dir}")
        else:
            file_dirs += glob.glob(os.path.join(dest_path, name), recursive=recursive)
            logger.info(f"find files: {len(file_dirs)}")
        
        file_nums = len(file_dirs)
        for file_filter in tqdm.tqdm(filter_files, total=len(filter_files), desc="Filter files"):
            if file_filter in file_dirs:
                file_dirs.remove(file_filter)

        print(f"sum file: {len(file_dirs)}; Filter files num: {file_nums - len(file_dirs)}")

        random.shuffle(file_dirs)
        train_data = file_dirs[:int(len(file_dirs) * train_ratio)]
        val_data = file_dirs[int(len(file_dirs) * train_ratio):int(len(file_dirs) * (train_ratio + val_ratio))]
        test_data = file_dirs[int(len(file_dirs) * (1 - test_ratio)):]
        # reference segment to split data of human seq
        print(f"train: val: test = {len(train_data)}: {len(val_data)}: {len(test_data)}")
        if os.path.exists(save_path) is False:
            os.makedirs(save_path, exist_ok=True)
            print(f"Create save path: {save_path}")
        self.write_txt(train_data, save_path, "train.txt")
        self.write_txt(val_data, save_path, "val.txt")
        self.write_txt(test_data, save_path, "test.txt")

    @staticmethod
    def _set_random_seed(seed):
        random.seed(seed)
        np.random.seed(seed)

    @staticmethod
    def write_txt(data: list, save_dir: str, name: str):
        if data:
            assert os.path.exists(save_dir), f"{save_dir} is not exit."
            path = os.path.join(save_dir, name)
            with open(path, "w", encoding="utf8") as f:
                f.write("\n".join(data) + "\n")
        else:
            print(f"ERROR: {name} is empty. data nums: {len(data)}")


def main(args):
    t2t_data = T2TDataProcess()
    t2t_data.split_pretrain_dta(
        dest_path=args.dest_path,
        save_path=args.save_path,
        train_ratio=args.train_ratio,
        val_ratio=args.val_ratio,
        test_ratio=args.test_ratio,
        name=args.file_name,
        recursive=args.recursive,
        filter_file_path=args.filter_file_path
    )


if __name__ == '__main__':
    args = parse_args()
    main(args)

    # dest_path = "/home/share/huadjyin/home/sunhaotong/trash/ll/cut_data/T2T/chm13v2.0_134144_no_lap"
    # save_path = "/home/share/huadjyin/home/s_sukui/02_data/07_genomics_data/T2T/pretrain"
    # t2t_data = T2TDataProcess()
    # t2t_data.split_pretrain_dta(dest_path, save_path)
