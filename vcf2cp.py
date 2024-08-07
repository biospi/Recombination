#!/usr/bin/env python3
import typer
import sys
import re
from tqdm import tqdm
from pathlib import Path


# Helper function to print help message and exit
def print_help():
    print("CONVERTS PHASED SHAPEIT/IMPUTE2 OUTPUT TO CHROMOPAINTER-STYLE INPUT FILES\n")
    print("usage:   python vcf2cp.py <options> input.vcf output_prefix\n")
    print("where:\n")
    print("        (i) input.vcf = phased vcf file with GT field.\n")
    print("        (ii) output_prefix = filename prefix for chromopainter input file(s). The suffix \".phase\" is added\n\n")
    print("The output is in CHROMOPAINTER v2 input format.\n")
    print("<options>:\n")
    print("-J:                 Jitter (add 1) snp locations if snps are not strictly ascending. By default, an error is produced.\n")
    print("-m <val>:           Maximum individuals processed in a single pass. Setting this larger is faster but uses more memory. Default: 500.\n")
    print("-g <val>:           Gap between chromosomes. By default, provides an error with multiple chromosome data, but you can set it to 100000000 or some other large value to make whole genome analysis possible.\n")
    print("-p <val>:           Ploidy. Default: 2.\n")
    print("NOTE: TO USE IN CHROMOPAINTER: You also need a recombination map. Create this with the \"convertrecfile.pl\" or \"makeuniformrecfile.pl\" scripts provided.\n\n")
    print(" !!! WARNING:  THIS PROGRAM DOES NOT SUFFICIENTLY CHECK FOR MISSPECIFIED FILES. WE ARE NOT ACCOUNTABLE FOR THIS RUNNING INCORRECTLY !!!\n")
    sys.exit()

# Function to process the VCF file
def process_vcf_header(vcf_file):
    ids = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#CHROM'):
                ids = line.strip().split('\t')[9:]
                break
    return ids

# Function to process the VCF file to identify SNPs
def process_vcf_snps(vcf_file, jitter):
    posvec = []
    nsnps = 0
    chrom = ""
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            items = line.strip().split('\t')
            nsnps += 1
            if chrom == "":
                chrom = items[0]
            else:
                if chrom != items[0] and chromosomegap < 0:
                    sys.exit("Require all SNPs to be on the same chromosome!")
            pos = int(items[1])
            if posvec and pos <= posvec[-1] and pos >= 0:
                if not jitter:
                    sys.exit("ERROR: SNPs are not strictly ascending, exiting. Rerun with -J to jitter the SNP locations, or remove multi-allelic SNPs.")
                print(f"Duplication found: setting {pos} to {posvec[-1] + 1}")
                pos = posvec[-1] + 1
            posvec.append(pos)
    return nsnps, posvec

# Function to process individuals and write output
def process_individuals(vcf_file, output_prefix, numindsmaxsize, nsnps, posvec, ploidy):
    if Path(f"{output_prefix}.phase").exists():
        print("output already exists.")
        return
    
    ids = process_vcf_header(vcf_file)
    ninds = len(ids)
    nhaps = ninds * ploidy
    numsplits = (ninds + numindsmaxsize - 1) // numindsmaxsize

    with open(f"{output_prefix}.ids", 'w') as out_ids:
        out_ids.write('\n'.join(ids) + '\n')

    with open(f"{output_prefix}.phase", 'w') as out_phase:
        out_phase.write(f"{nhaps}\n")
        out_phase.write(f"{nsnps}\n")
        out_phase.write("P " + ' '.join(map(str, posvec)) + '\n')

        for a in tqdm(range(numsplits)):
            start_ind = a * numindsmaxsize
            end_ind = min((a + 1) * numindsmaxsize, ninds)
            genomat = [[None] * nsnps for _ in range(ploidy * (end_ind - start_ind))]
            with open(vcf_file, 'r') as f:
                snpon = 0
                for line in f:
                    if line.startswith('#'):
                        continue
                    items = line.strip().split('\t')
                    for n in range(start_ind, end_ind):
                        gt = items[9 + n].split(':')[0]
                        gtv = gt.split('|')
                        if len(gtv) != ploidy:
                            sys.exit(f"Individual {ids[n]} with field {gt} at position index {snpon} not phased or diploid?")
                        for p in range(ploidy):
                            genomat[(n - start_ind) * ploidy + p][snpon] = gtv[p]
                    snpon += 1

            for hap in genomat:
                out_phase.write(''.join(hap) + '\n')

def vcf2cp(
    input_vcf: str = typer.Argument(..., help="Phased VCF file with GT field."),
    output_prefix: str = typer.Argument(..., help="Filename prefix for ChromoPainter input files."),
    jitter: bool = typer.Option(False, "--J", "-J", help="Jitter (add 1) SNP locations if SNPs are not strictly ascending."),
    numindsmaxsize: int = typer.Option(500, "--m", "-m", help="Maximum individuals processed in a single pass. Default: 500."),
    chromosomegap: int = typer.Option(-1, "--g", "-g", help="Gap between chromosomes. Default: -1."),
    ploidy: int = typer.Option(2, "--p", "-p", help="Ploidy. Default: 2.")
):
    output_prefix = re.sub(r'\.phase$', '', output_prefix)

    print("Options in effect:")
    print(f"  input file: {input_vcf}")
    print(f"  output ids file: {output_prefix}.ids")
    print(f"  output phase file: {output_prefix}.phase")
    print(f"  ploidy: {ploidy}")
    print(f"  jitter: {jitter}")
    print(f"  individuals per pass: {numindsmaxsize}")
    print(f"  chromosome gap (-ve for single chromosome analysis only): {chromosomegap}")

    #ids = process_vcf_header(input_vcf)
    nsnps, posvec = process_vcf_snps(input_vcf, jitter)
    process_individuals(input_vcf, output_prefix, numindsmaxsize, nsnps, posvec, ploidy)

if __name__ == "__main__":
    typer.run(vcf2cp)