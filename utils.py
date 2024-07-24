
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



def get_header(file_path, thresh=10):
    output_file_path = Path("header.txt")
    if output_file_path.exists():
        output_file_path.unlink()

    lines = []
    with gzip.open(file_path, 'rt') as f:
        for i, line in enumerate(f):
            #print(f"Line {i + 1} length: {len(line.strip())}")
            line = line.replace('\x00', '')
            if len(line.strip()) > 1500:
                break
            line = line.replace('\t',',').strip()
            #print(line)
            if 'na' in line.lower():
                continue
            if '(' in line or ')' in line:
                continue
            if len(line) < 1:
                continue
            lines.append(line)
    

    #TODO fix header parsing not generalisable!
    cleaned_header = lines[0].replace('\x00', '')
    cleaned_header = cleaned_header.replace('%', '')
    cleaned_header = f"sample_id{cleaned_header.split('sample_id')[-1]}"
    print(f"cleaned_header: {cleaned_header}")
    lines[0] = cleaned_header
        
    print(f"number of lines in header: {len(lines)}")
    with open(output_file_path.as_posix(), 'w') as file:
        for line in lines:
            file.write(line + '\n')

    cols = cleaned_header.split(',')
    print(f"cols: {cols}")
    df = pd.DataFrame([x.split(',') for x in lines[1:]], columns=cols)
    print(df)
    df.to_csv("header.csv", index=False)
    print(output_file_path)

    df_ = df[['sample_id', 'intergenic', 'indels', 'insertion', 'deletion', 'intergenic_indels','genic_indels']]
    values = df_.drop(columns=['sample_id'])
    print(f"value:{values}")
    values = values.astype(float)
    mask = (values > thresh).all(axis=1)
    result = df_.loc[mask, 'sample_id']
    id_to_exluce = result.values.tolist()
    df_id = pd.DataFrame([[x] for x in id_to_exluce])
    print(df_id)
    df_id.to_csv("id_to_exclude.csv", index=False)
    print(f"n id_to_exluce:{len(id_to_exluce)}")
    return id_to_exluce


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
                            if str(record.id).strip().lower() in [str(x).strip().lower for x in exlude]:
                                print(f"Excluding sequence {record.id} with length {seq_len}")
                                continue

                            seq_len = len(record.seq)
                            char_count = Counter(record.seq)  # Count each character in the sequence
                            log = f"{seq_len},{expected_length},{record.id},{record},{str(dict(char_count)).replace(',', ' ')}"
                            # print(log)
                            logs.append(log)
                            sequences.append(record)

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


                        if len(sequences) > 0:
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


if __name__ == "__main__":
    test_string = "wgs-mapping/0000775000175000017500000000000014502606564012747 5ustar  ubuntuubuntuwgs-mapping/mapping_summary-NC_011294.tsv0000664000175000017500000512065314020412051020024 0ustar  ubuntuubuntusample_id,filesize_R1,filesize_R2,reads,reference,ref_length,mapped_reads,%mapped,unmapped_reads,%unmapped,average_coverage,Ns,%Ns,variants,snps,nonsyn,syn,intergenic,indels,insertion,deletion,intergenic_indels,genic_indels"
    test_string = test_string.replace("%", '')
    test_string = f"sample_id{test_string.split('sample_id')[-1]}"
    print(test_string)