
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
import ast

import gzip
from mimetypes import guess_type
from functools import partial
from Bio import SeqIO

from tqdm import tqdm

def get_header(output_dir, file_path, thresh=10):
    output_file_path = output_dir / "header.txt"
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
    df.to_csv(output_dir / "header.csv", index=False)
    print(output_file_path)
    # df["mapped"] = df["mapped"].astype(float)
    # print(f"mapped:{df['mapped']}")
    # df_to_exlude = df[df["mapped"] < 98][["sample_id", "mapped"]]
    # print(df_to_exlude)
    # df_to_exlude.to_csv("id_to_exclude.csv", index=False)
    # print(f"n id_to_exluce:{len(df_to_exlude)}/{len(df['sample_id'].values.tolist())}")
    # ids_to_exlude = df_to_exlude["sample_id"].values.tolist()
    # #print(ids_to_exlude)

    # df_to_keep = df[df["mapped"] > 98][["sample_id", "mapped"]]
    # df_to_keep.to_csv("id_to_keep.csv", index=False)
    # ids_to_keep = df_to_keep["sample_id"].values.tolist()
    return 


def get_consec(char,record, thresh=1000):
    # Find all occurrences of consecutive char
    consecutive_Ns = re.finditer(f'{char}+', str(record.seq))
    # Store each occurrence and its length
    counts = np.array([len(match.group()) for match in consecutive_Ns])
    counts = counts[counts > thresh]
    return counts
    

def filter_isolates_in_tar_gz(out_dir, tar_gz_filepath, cleaned_filepath, to_keep):
    print(f"Cleaning raw data in {cleaned_filepath}...")
    filter = [str(x).strip().lower().replace('"', '').replace("'", '') for x in to_keep]
    #print(f"filter={filter}")
    #sequence_lengths = []
    expected_length = -1
    logs = []
    
    # First, collect all sequence lengths
    # with tarfile.open(tar_gz_filepath, "r:gz") as tar:
    #     for member in tar.getmembers():
    #         if member.isfile() and member.name.endswith('.fasta'):
    #             file = tar.extractfile(member)
    #             if file is not None:
    #                 with TextIOWrapper(file, encoding='utf-8') as fasta_file:
    #                     for record in SeqIO.parse(fasta_file, "fasta"):
    #                         sequence_lengths.append(len(record.seq))
    
    # # Determine the expected length as the most common sequence length
    # df = pd.DataFrame(sequence_lengths, columns=['Length'])
    # expected_length = df['Length'].mode()[0]
    # print(f"Expected sequence length: {expected_length}")

    # Now filter sequences based on the expected length and minimum length
    with tarfile.open(tar_gz_filepath, "r:gz") as tar, open(cleaned_filepath, "w") as cleaned:
        for member in tar.getmembers():
            if member.isfile() and member.name.endswith('.fasta'):
                file = tar.extractfile(member)
                if file is not None:
                    with TextIOWrapper(file, encoding='utf-8') as fasta_file:
                        sequences = []
                        cpt = 0
                        for record in SeqIO.parse(fasta_file, "fasta"):
                            r = str(record.id).strip().lower().replace('"', '').replace("'", '')
                            if r not in filter:
                                #print(f"Excluding sequence {record.id}")
                                continue

                            seq_len = len(record.seq)
                            char_count = Counter(record.seq)  # Count each character in the sequence
                            #print(char_count)
                            if char_count["N"] > 3500:
                                continue

                            log = f"{seq_len},{expected_length},{record.id},{record},{str(dict(char_count)).replace(',', ' ')}"
                            # print(log)
                            logs.append(log)
                            sequences.append(record)
                            cpt += 1
                            SeqIO.write([record], cleaned, "fasta")

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


                        # if len(sequences) > 0:
                        #     print(f"Found {len(sequences)} valid sequences. Exporting to new archive...")
                        #     # Create a temporary file to write filtered sequences
                        #     tmp_file = StringIO()
                        #     SeqIO.write(sequences, tmp_file, "fasta")
                        #     tmp_file.seek(0)

                        #     # Convert StringIO content to bytes
                        #     byte_file = BytesIO(tmp_file.getvalue().encode('utf-8'))
                        #     byte_file.seek(0)

                        #     # Create a new tarinfo object for the filtered sequences file
                        #     tarinfo = tarfile.TarInfo(name=member.name)
                        #     tarinfo.size = len(byte_file.getvalue())
                        #     byte_file.seek(0)  # Reset the pointer to the beginning

                        #     # Add the filtered sequences file to the new tar.gz archive
                        #     cleaned_tar.addfile(tarinfo, fileobj=byte_file)

    # Save logs to a CSV file
    log_df = pd.DataFrame([log.split(",") for log in logs], columns=["seq_len", "expected_length", "record_id", "record", "count"])
    log_csv_filepath = out_dir / "log.csv"
    log_df.to_csv(log_csv_filepath, index=False)
    print(f"Logs saved to {log_csv_filepath}")
    return cpt


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


def run_cmd(cmd_str, out_dir, tag, n_isolate=-1):
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


def convert_to_dict(count_str):
    count_str = count_str.replace("  ", ", ")
    return ast.literal_eval(count_str)


def add_missigness(input_data):
    print(input_data)
    df = pd.read_csv(input_data)
    df["missingness__autocolour"] = np.nan
    print(df.columns)

    df_log = pd.read_csv(Path('output/log.csv'))
    df_log['count'] = df_log['count'].apply(convert_to_dict)
    df_log['missingness'] = df_log['count'].apply(lambda x: dict(x).get('N', 0))
    df_log = df_log[['record_id', 'missingness']]
    print(df_log)

    for _, row in df_log.iterrows():
        df.loc[df["id"] == row['record_id'], "missingness__autocolour"] = row["missingness"]
    
    df.to_csv("data_with_miss.csv", index=False)


def keep_top(file: Path, output_file: Path):
    with file.open('r') as infile:
        with output_file.open('w') as outfile:
            count = 0
            for line in infile:
                outfile.write(line)
                count += 1
                if count >= 1000:
                    break


def append_fasta(fasta_file, data_to_append, output_filepath):
    print(f"Append fasta {fasta_file}...")
    logs = []
    with open(output_filepath, "w") as out_file:
        cpt = 0
        if cpt == 0:
            #read data_to_append file and write                 
            for record in SeqIO.parse(data_to_append, "fasta"):
                record.description = record.description.replace("\n",".").replace(",", "")
                seq_len = len(record.seq)
                char_count = Counter(record.seq)
                SeqIO.write([record], out_file, 'fasta')
                r_f = str(record).replace("\n",".").replace(",", "")
                log = f"{seq_len},{record.id},{r_f},{str(dict(char_count)).replace(',', ' ')}"
                print(log)
                logs.append(log)
                cpt+=1

        for record in SeqIO.parse(fasta_file, "fasta"):
            record.description = record.description.replace("\n",".").replace(",", "")
            seq_len = len(record.seq)
            char_count = Counter(record.seq)  # Count each character in the sequence
            #print(char_count)
            if char_count["N"] > 80000:
                print(f"remove {record.id}....")
                continue

            log = f"{seq_len},{record.id},{record},{str(dict(char_count)).replace(',', ' ')}"
            # print(log)
            logs.append(log)
            cpt += 1
            SeqIO.write([record], out_file, 'fasta')
            # if cpt > 100:
            #     break
    print(f"Number of isolates left {cpt}.")
    # Save logs to a CSV file
    log_df = pd.DataFrame([log.split(",") for log in logs], columns=["seq_len", "record_id", "record", "count"])
    log_csv_filepath = data_to_append.parent / "log_4.csv"
    log_df.to_csv(log_csv_filepath, index=False)
    print(f"Logs saved to {log_csv_filepath}")
    return cpt


def calculate_missingness(group, sample_start_index):
    total_samples = 0
    missing_samples = 0

    for line in group:
        columns = line.strip().split("\t")
        # Process each sample in the VCF (starting from the 9th column onwards)
        for sample_data in columns[sample_start_index:]:
            total_samples += 1
            # Check if the sample has missing data ('*' or 'N')
            if sample_data == '*' or sample_data == 'N':
                missing_samples += 1

    # Calculate missingness percentage
    missingness_percentage = (missing_samples / total_samples) * 100
    return missingness_percentage


# def filter_vcf(input_vcf, output_vcf, missingness_threshold=10.0):
#     print(f"filter vcf missingness_threshold={missingness_threshold}")
#     current_pos = None
#     group = []
#     num_lines = sum(1 for line in open(input_vcf))

#     with open(input_vcf, 'r') as vcf_in, open(output_vcf, 'w') as vcf_out:
#         for line in tqdm(vcf_in, total=num_lines):
#             if line.startswith("#"):
#                 vcf_out.write(line)
#                 # Save the index from where the sample data starts
#                 if line.startswith("#CHROM"):
#                     headers = line.strip().split("\t")
#                     sample_start_index = 9  # Sample data starts from the 9th column
#             else:
#                 columns = line.strip().split("\t")
#                 pos = columns[1]  # The position of the current line

#                 if current_pos is None:
#                     current_pos = pos

#                 if pos == current_pos:
#                     group.append(line)
#                 else:
#                     missingness_percentage = calculate_missingness(group, sample_start_index)
#                     if missingness_percentage <= missingness_threshold:
#                         for group_line in group:
#                             vcf_out.write(group_line)
#                     else:
#                         print(f"Remove pos={pos}")
#                     group = [line]
#                     current_pos = pos
#         if group:
#             missingness_percentage = calculate_missingness(group, sample_start_index)
#             if missingness_percentage <= missingness_threshold:
#                 for group_line in group:
#                     vcf_out.write(group_line)
#     print(output_vcf)


def filter_samples_by_missingness_and_remove_star(vcf_file, output_file, threshold=0.1):
    with open(vcf_file, 'r') as f:
        lines = f.readlines()
    
    header = []
    genotype_data = []
    num_variants = 0
    for line in lines:
        if line.startswith("#"):
            header.append(line)
        else:
            row = line.strip().split('\t')
            if '*' not in row[4]: 
                num_variants += 1
                genotype_data.append(row)

    sample_columns = [list(x) for x in zip(*genotype_data)][9:]  
    
    total_variants = len(genotype_data)
    samples_to_keep = []
    
    for i, sample in enumerate(sample_columns):
        missing_genotypes = sum(1 for gt in sample if gt == '.')
        missingness = missing_genotypes / total_variants
        
        if missingness <= threshold:
            samples_to_keep.append(i + 9)  
    
    with open(output_file, 'w') as out:
        # Write the header
        out.write(''.join(header[:9])) 
        out.write('\t'.join([header[9].split()[i] for i in samples_to_keep]) + '\n') 

        for line in genotype_data:
            out.write('\t'.join(line[:9]))  
            out.write('\t'.join([line[i] for i in samples_to_keep]) + '\n')  


def build_target_files(vcf_file):
    vcf_file = Path(vcf_file)
    out_dir = vcf_file.parent / "targets"
    out_dir.mkdir(parents=True, exist_ok=True)
    total_lines = sum(1 for _ in open(vcf_file))

    with vcf_file.open() as f:
        header = []
        sample_names = []
        sample_columns = {}
        variant_counts = {}  # To store variability levels

        # Parse the VCF file line by line
        for line in tqdm(f, total=total_lines):
            if line.startswith("##"):  # VCF metadata header lines
                header.append(line)
            elif line.startswith("#CHROM"):  # The column names (including samples)
                header.append(line)
                columns = line.strip().split("\t")
                sample_names = columns[9:]  # Sample names start at column 9
                sample_columns = {name: [] for name in sample_names}
                variant_counts = {name: {"total": 0, "variant": 0} for name in sample_names}
            else:
                # Extract the variant information and corresponding sample columns
                variant_info = line.strip().split("\t")
                for i, sample_name in enumerate(sample_names):
                    sample_data = variant_info[9 + i]
                    # Add the variant info plus the sample data
                    sample_columns[sample_name].append("\t".join(variant_info[:9] + [sample_data]))
                    
                    # Calculate variability level
                    variant_counts[sample_name]["total"] += 1
                    if sample_data.startswith("1") or sample_data.startswith("0/1"):
                        variant_counts[sample_name]["variant"] += 1

    # Create individual VCF files for each sample
    for sample_name in sample_names:
        total_variants = variant_counts[sample_name]["total"]
        non_ref_variants = variant_counts[sample_name]["variant"]
        # Calculate variability level as a percentage
        variability_level = int((non_ref_variants / total_variants) * 100) if total_variants > 0 else 0
        
        # Prepend variability level to the filename
        target_file = out_dir / f"{variability_level}_{sample_name}.vcf"
        print(target_file)
        with target_file.open("w") as f:
            # Write the header and then the variant information for the sample
            f.writelines(header)
            f.write("\n".join(sample_columns[sample_name]) + "\n")

    print(f"Target files created in {out_dir}")


if __name__ == "__main__":
    # i = 0
    # for line in open(Path('output/wgs-mapping_with_ref.fasta.expand.vcf')):
    #     print(line[0:50])
    #     i += 1
    #     if i > 20:
    #         break

    # test_string = "wgs-mapping/0000775000175000017500000000000014502606564012747 5ustar  ubuntuubuntuwgs-mapping/mapping_summary-NC_011294.tsv0000664000175000017500000512065314020412051020024 0ustar  ubuntuubuntusample_id,filesize_R1,filesize_R2,reads,reference,ref_length,mapped_reads,%mapped,unmapped_reads,%unmapped,average_coverage,Ns,%Ns,variants,snps,nonsyn,syn,intergenic,indels,insertion,deletion,intergenic_indels,genic_indels"
    # test_string = test_string.replace("%", '')
    # test_string = f"sample_id{test_string.split('sample_id')[-1]}"
    # print(test_string)
    #add_missigness(Path('output/data.csv'))
    #keep_top(Path('output/wgs-mapping_with_ref.fasta.vcf'), Path('top_1000_wgs-mapping_with_ref.fasta.vcf'))
    #filter_samples_by_missingness(Path('output/wgs-mapping_with_ref.fasta.vcf'), Path('output/wgs-mapping_with_ref.filtered.fasta.vcf'))
    build_target_files(Path('output/wgs-mapping.tar.gz.filtered.ref.vcf'))

    # append_fasta(Path("wgs-mapping/wgs-mapping/mapping-NC_011294.filtered.fa"), 
    #              Path("output/S_enterica_Enteritidis_P12510.fasta"), 
    #              Path("output/wgs-mapping_with_ref.fasta"))