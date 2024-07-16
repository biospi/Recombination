
import subprocess
import time
import datetime
from Bio import SeqIO
import pandas as pd
import tarfile
import gzip
from io import TextIOWrapper
from Bio import SeqIO
from pathlib import Path


def count_isolates_in_tar_gz(tar_gz_filepath):
    count = 0
    
    with tarfile.open(tar_gz_filepath, "r:gz") as tar:
        for member in tar.getmembers():
            if member.isfile() and member.name.endswith('.fasta'):
                file = tar.extractfile(member)
                if file is not None:
                    with TextIOWrapper(file, encoding='utf-8') as fasta_file:
                        for record in SeqIO.parse(fasta_file, "fasta"):
                            count += 1
    return count


def clean_data(input_file):
    print("cleaning data...")
    df = pd.read_csv(input_file, nrows=10, on_bad_lines='warn')
    print(df)
    return df


def count_isolates(fasta_filepath):
    count = 0
    with open(fasta_filepath, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            count += 1
    print(f"Found {count} isolates in {fasta_filepath}.")
    return count


def run_cmd(cmd_str, out_dir, tag, n_isolate):
    print(cmd_str)
    start = time.time()
    result = subprocess.call(cmd_str, shell=True, stdout=subprocess.PIPE)
    print(f"cmd status:{result}")
    dt = datetime.timedelta(seconds=time.time() - start)
    print(f"total time (H:M:S.mm)= {str(dt)}")
    filedir = out_dir / "time"
    filedir.mkdir(parents=True, exist_ok=True)
    filename = f"{tag}.csv"
    filepath = filedir / filename
    df = pd.DataFrame()
    df["cmd_str"] = [cmd_str]
    df["time_s"] = [dt.total_seconds()]
    df["time_m"] = [dt.microseconds]
    df["n_isolate"] = [n_isolate]
    df.to_csv(filepath, columns=['cmd_str', 'time_s', 'time_m', 'n_isolate'], index=False)
    return cmd_str