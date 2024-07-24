
import subprocess
import time
import datetime
from Bio import SeqIO
import pandas as pd
import tarfile
import gzip
from io import TextIOWrapper, BytesIO, StringIO
from Bio import SeqIO
from pathlib import Path
from collections import Counter
import re
import numpy as np
import gzip


def get_header(file_path):
    output_file_path = Path("header.csv")
    if output_file_path.exists():
        output_file_path.unlink()

    lines = []
    with gzip.open(file_path, 'rt') as f:
        for i, line in enumerate(f):
            #print(f"Line {i + 1} length: {len(line.strip())}")
            if len(line.strip()) > 1500:
                break
            line = line.replace('\t',',').strip()
            print(line)
            lines.append(line)
    print(f"number of lines in header: {len(lines)}")
    df = pd.DataFrame(lines, columns=lines[0])
    df.to_csv(output_file_path, index=False)
    print(output_file_path)


def get_consec(char,record, thresh=1000):
    # Find all occurrences of consecutive char
    consecutive_Ns = re.finditer(f'{char}+', str(record.seq))
    # Store each occurrence and its length
    counts = np.array([len(match.group()) for match in consecutive_Ns])
    counts = counts[counts > thresh]
    return counts
    

def filter_isolates_in_tar_gz(tar_gz_filepath, cleaned_tar_gz_filepath, exlude):
    print(f"Cleaning raw data in {cleaned_tar_gz_filepath}...")
    print(f"exlude={exlude}")
    sequence_lengths = []
    logs = []
    
    # First, collect all sequence lengths
    with tarfile.open(tar_gz_filepath, "r:gz") as tar:
        for member in tar.getmembers():
            if member.isfile() and member.name.endswith('.fasta'):
                file = tar.extractfile(member)
                if file is not None:
                    with TextIOWrapper(file, encoding='utf-8') as fasta_file:
                        for record in SeqIO.parse(fasta_file, "fasta"):
                            sequence_lengths.append(len(record.seq))
    
    # Determine the expected length as the most common sequence length
    df = pd.DataFrame(sequence_lengths, columns=['Length'])
    expected_length = df['Length'].mode()[0]
    print(f"Expected sequence length: {expected_length}")

    # Now filter sequences based on the expected length and minimum length
    with tarfile.open(tar_gz_filepath, "r:gz") as tar, tarfile.open(cleaned_tar_gz_filepath, "w:gz") as cleaned_tar:
        for member in tar.getmembers():
            if member.isfile() and member.name.endswith('.fasta'):
                file = tar.extractfile(member)
                if file is not None:
                    with TextIOWrapper(file, encoding='utf-8') as fasta_file:
                        sequences = []
                        for record in SeqIO.parse(fasta_file, "fasta"):
                            #print(record.id)
                            # if record.id not in ['SRR5193283', 'SRR3049562', 'SRR6900352']:
                            #     continue


                            seq_len = len(record.seq)
                            char_count = Counter(record.seq)  # Count each character in the sequence
                            log = f"{seq_len},{expected_length},{record.id},{record},{str(dict(char_count)).replace(',', " ")}"
                            print(log)
                            logs.append(log)

                            # A_counts = get_consec('A', record)
                            # G_counts = get_consec('G', record)
                            # T_counts = get_consec('T', record)
                            # N_counts = get_consec('N', record)

                            # if record.id in exlude or len(A_counts) > 0 or len(G_counts) > 0 or len(T_counts) > 0 or len(N_counts) > 0:
                            #     print(f"Excluding sequence {record.id} with length {seq_len} A_counts={A_counts}, G_counts={G_counts}, T_counts={T_counts}, N_counts={N_counts}")
                            #     continue

                            # if seq_len == expected_length:
                            #     sequences.append(record)
                            # else:
                            #     print(f"Excluding sequence {record.id} with length {seq_len}")

                        if sequences:
                            # Create a temporary file to write filtered sequences
                            tmp_file = StringIO()
                            SeqIO.write(sequences, tmp_file, "fasta")
                            tmp_file.seek(0)

                            # Convert StringIO content to bytes
                            byte_file = BytesIO(tmp_file.getvalue().encode('utf-8'))
                            byte_file.seek(0)

                            # Create a new tarinfo object for the filtered sequences file
                            tarinfo = tarfile.TarInfo(name=member.name)
                            tarinfo.size = len(byte_file.getvalue())
                            byte_file.seek(0)  # Reset the pointer to the beginning

                            # Add the filtered sequences file to the new tar.gz archive
                            cleaned_tar.addfile(tarinfo, fileobj=byte_file)

    # Create a DataFrame from the lengths list
    length_counts = df['Length'].value_counts().reset_index()
    length_counts.columns = ['Length', 'Count']
    lengths_csv_filepath = Path("sequence_lengths.csv")
    length_counts.to_csv(lengths_csv_filepath, index=False)
    print(f"Sequence lengths saved to {lengths_csv_filepath}")

    # Save logs to a CSV file
    log_df = pd.DataFrame([log.split(",") for log in logs], columns=["seq_len", "expected_length", "record_id", "record", "count"])
    log_csv_filepath = Path("log.csv")
    log_df.to_csv(log_csv_filepath, index=False)
    print(f"Logs saved to {log_csv_filepath}")


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