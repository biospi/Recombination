import pandas as pd
import numpy as np

from makeuniformrecfile import makeuniformrecfile
from utils import run_cmd
from vcf2cp import vcf2cp
from visu import histogram_vcf
np.random.seed(0)
from pathlib import Path
from tqdm import tqdm

from collections import Counter


def clean_vcf(vcf_path, thresh=5):
    print(vcf_path)
    line_count = sum(1 for line in open(vcf_path))
    clean_vcf = vcf_path.parent / f"{vcf_path.stem}.filtered.vcf"
    if clean_vcf.exists():
        clean_vcf.unlink()
    print(clean_vcf)
    print("cleaning vcf file...")
    with open(vcf_path, 'r') as file:
        with open(clean_vcf, "w") as out_file:
            for i, line in tqdm(enumerate(file), total=line_count):
                if line.startswith('#CHROM'): 
                    columns = line.strip().split('\t')
                    n_samples = len(columns) - 9 
                if i >= 4:
                    split = line.split('\t')
                    new_line = '\t'.join(split)
                    if '*' in new_line:
                        #print('\t'.join([x[0:4] for x in new_line.split('\t')][0:20]))  
                        #char_count = Counter(new_line)
                        genotypes = split[9:] 
                        char_count = Counter(genotypes)
                        missing_count = char_count['1']
                        missing_percentage = (missing_count / n_samples) * 100

                        if missing_percentage > thresh:
                            # print('\t'.join([x[0:4] for x in new_line.split('\t')][0:20]))  
                            # print(char_count)
                            # print("MISSING!!!!")
                            continue
                        out_file.write(new_line)
                    else:
                        out_file.write(line)
                else:
                    out_file.write(line)
    print(f"clean_vcf={clean_vcf}")
    return clean_vcf


def build_ref_from_vcf(vcf_path):
    print(vcf_path)
    line_count = sum(1 for line in open(vcf_path))
    ref_vcf = vcf_path.parent / f"{vcf_path.stem}.ref.vcf"
    if ref_vcf.exists():
        ref_vcf.unlink()
    print(ref_vcf)
    print("build ref vcf file...")
    with open(vcf_path, 'r') as file:
        with open(ref_vcf, "w") as out_file:
            for i, line in tqdm(enumerate(file), total=line_count):
                #print('\t'.join([x[:] for x in line.split('\t')][0:20]))  
                # if i > 100:
                #     break
                if i >= 4:
                    split = line.split('\t')
                    # split[np.random.randint(9, 20)] = '2'
                    # split[np.random.randint(9, 20)] = '1'
                    new_line = '\t'.join(split)
                    if '*' in new_line:
                        #print('\t'.join([x[0:4] for x in new_line.split('\t')][0:20]))  
                        if split[4][0] == '*':
                            new_line = '\t'.join(new_line.split('\t')[0:9]) + '\t' + '\t'.join(new_line.split('\t')[9:]).replace('1', '0').replace('2', '1')
                        if split[4][-1] == '*':
                            new_line = '\t'.join(new_line.split('\t')[0:9]) + '\t' + '\t'.join(new_line.split('\t')[9:]).replace('2', '0')
                        new_line = new_line.replace('*,', '')
                        new_line = new_line.replace(',*', '')
                        #new_line = new_line.replace('*', '.')
                        #print('\t'.join([x[0:4] for x in new_line.split('\t')][0:20]))  
                        out_file.write(new_line)
                    else:
                        #print('\t'.join([x[0:4] for x in new_line.split('\t')][0:20]))  
                        out_file.write(line)
                else:
                    out_file.write(line)
    print(f"ref_vcf={ref_vcf}")
    return ref_vcf



if __name__ == "__main__":
    #filepath_vcf_expand = Path('output/wgs-mapping_with_ref.fasta.vcf')
    # output_dir = Path("output2") / 'ref'
    # output_dir.mkdir(parents=True, exist_ok=True)
    # build_ref_from_vcf(filepath_vcf_expand)

    vcf_filepath = Path('output/wgs-mapping.tar.gz.vcf')
    vcf_filepath_ref = Path('output/wgs-mapping.tar.gz.filtered.ref.vcf')
    # histogram_vcf(vcf_filepath, title="(Raw)")
    # filtered_vcf_filepath = clean_vcf(vcf_filepath, thresh=10)
    # histogram_vcf(filtered_vcf_filepath, title="(Filtered)")
    # vcf_filepath_ref = build_ref_from_vcf(filtered_vcf_filepath)

    sites_filepath = vcf_filepath.parent / "sites.txt"
    map_filepath = vcf_filepath.parent / "map.txt"
    # run_cmd(
    #     f"bcftools query -f '%POS\\n' {vcf_filepath} > {sites_filepath}",
    #     vcf_filepath.parent,
    #     "bcftools",
    #     -1,
    # )
    # run_cmd(
    #     f"echo -e 'pd\\tgd' > {map_filepath}",
    #     vcf_filepath.parent,
    #     "echo",
    #     -1,
    # )
    # run_cmd(
    #     f"awk '{{print $1\"\\t\"1e-6*$1}}' {sites_filepath.as_posix()} >> {map_filepath.as_posix()}",
    #     vcf_filepath.parent,
    #     "awk",
    #     -1,
    # )


    run_cmd(
        "cd /home/axel/tools/SparsePainter",
        vcf_filepath.parent,
        "cd",
        -1,
    )

    data_dir = "/home/axel/python-projects/Recombination/"

    run_cmd(
        f"./SparsePainter -reffile {data_dir}{vcf_filepath_ref.as_posix()} -targetfile {data_dir}{vcf_filepath_ref.as_posix()} -mapfile {data_dir}{map_filepath.as_posix()} -popfile {} -prob -haploid -chunklength -probstore raw",
        vcf_filepath.parent,
        "cd",
        -1,
    )


    # vcf2cp(input_vcf=filepath_vcf_expand.as_posix(), 
    #     output_prefix=(output_dir / "vcf2cp_").as_posix(),
    #     jitter=True,
    #     numindsmaxsize=500,
    #     chromosomegap=-1,
    #     ploidy=1)
    
    # makeuniformrecfile(phasefile=(output_dir / "vcf2cp_.phase").as_posix(),
    #                 outputfile=(output_dir / "vcf2cp_.recomb").as_posix())