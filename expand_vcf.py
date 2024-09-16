#!/usr/bin/env python3

import typer
import sys


def expand(
    input: str = typer.Option(..., "-i", "--input", help="Input PIRATE.*.tsv file [required]"),
    output: str = typer.Option(..., "-o", "--output", help="Output treeWAS input file [required]"),
    freq_low: float = typer.Option(0.00, "--fl", "--freq-low", help="Min SNP frequency to include in output [default: 0.00]"),
    freq_high: float = typer.Option(1.00, "--fh", "--freq-high", help="Max SNP frequency to include in output [default: 1.00]"),
    n_sites: bool = typer.Option(False, "-n", "--n-sites", help="Include sites called due to the presence of gaps [default: exclude]"),
    store_inv: bool = typer.Option(False, "-s", "--store-inv", help="Store alleles that are the inverse of another [default: off]"),
    sample_list: str = typer.Option(None, "-l", "--list", help="List of samples to include [default: off]"),
    freq_inv: float = typer.Option(0.98, "--freq-inv", help="Frequency threshold to be considered inverse allele [default: 0.98]"),
    data_type: str = typer.Option("auto", "-t", "--type", help="Type of data [default: auto, opts: nuc, aa]"),
):
    # Validations
    if not 0 <= freq_inv <= 1:
        print(" - ERROR: inverse allele frequency threshold (--freq-inv) should be a number between 0-1.")
        return -1
    if data_type not in ["nuc", "aa", "auto"]:
        print(" - ERROR: --type only accepts nuc, aa, auto.")
        return -1

    # Feedback - settings
    print("\nSettings:")
    print(" - including gaps in output." if n_sites else " - excluding gaps in output.")
    print(" - including inverse_alleles in output." if store_inv else " - excluding inverse_alleles in output.")
    print(f" - inverse allele frequency = {freq_inv}")

    # Test type of input
    aa = 0
    if data_type == "auto":
        test_count = 0
        non_atcg = 0
        try:
            with open(input, "r") as vcf:
                for line in vcf:
                    if line.startswith("#CHROM\tPOS") or line.startswith("#"):
                        continue
                    test_count += 1
                    ref = line.split("\t")[3].upper()
                    alt = line.split("\t")[4].upper()
                    alleles = [ref] + alt.split(",")
                    non_atcg = any(not set(allele).issubset({"A", "T", "C", "G", "N", "*"}) for allele in alleles)
                    if non_atcg:
                        aa = 1
                        break
                    if test_count >= 1000:
                        break
        except Exception as e:
            print(f"VCF file did not open: {e}")
            raise typer.Exit(code=1)

    print(" - vcf being processed as including amino acid sequence" if aa else " - vcf being processed as including nucleotide sequence")
    print("\nRunning:")

    # Open list file if provided
    sample_dict = {}
    no_samples_list = 0
    if sample_list:
        try:
            with open(sample_list, "r") as list_file:
                for line in list_file:
                    sample_id = line.strip().split("\t")[0]
                    if sample_id != "id":
                        sample_dict[sample_id] = 1
        except Exception as e:
            print(f" - ERROR: could not open list ({sample_list}): {e}")
            raise typer.Exit(code=1)
        no_samples_list = len(sample_dict)
        print(f" - {no_samples_list} samples to include from list file.")

    # Parse VCF file
    no_sites = 0
    no_variants = 0
    no_stored_vars = 0
    n_vars = 0
    headers = []
    include = []
    inverse = 0

    try:
        with open(input, "r") as vcf, open(output, "w") as out_file:
            for line in vcf:
                if line.startswith("#CHROM\tPOS"):
                    headers = line.strip().split("\t")
                    if sample_list:
                        include = [i for i, h in enumerate(headers) if h in sample_dict]
                        for sample in sample_dict:
                            if sample_dict[sample] == 1:
                                print(f"missing samples:\n{sample}")
                    else:
                        include = list(range(9, len(headers)))
                    out_file.write("\t".join(headers[:9] + [headers[i] for i in include]) + "\n")
                    no_samples = len(include)
                    if sample_list:
                        print(f" - {no_samples} samples of {no_samples_list} found in vcf headers.")
                    else:
                        print(f" - {no_samples} samples found in vcf headers.")
                    print(" - 0 variant sites processed.")
                elif line.startswith("#"):
                    out_file.write(line)
                else:
                    no_sites += 1
                    line_split = line.strip().split("\t")
                    pos = line_split[1]
                    ref = line_split[3].upper().replace("-", "*").replace("N", "*" if aa == 0 else "N")
                    alt = line_split[4].upper().replace("-", "*").replace("N", "*" if aa == 0 else "N")
                    alleles = [ref] + alt.split(",")
                    n_Ns = alleles.count("*")
                    n_vars += n_Ns
                    no_alleles = len(alleles) - (0 if n_sites else n_Ns)
                    no_variants += len(alleles)

                    site_store = {}
                    count_store = {}

                    if no_alleles > 1:
                        for a_idx, a_char in enumerate(alleles):
                            if n_sites or (not n_sites and a_char != "*"):
                                outline = line_split[:4] + [a_char] + line_split[5:9]
                                s_count = sum(1 for i in include if line_split[i] == str(a_idx))
                                freq = s_count / no_samples
                                if freq_low <= freq <= freq_high:
                                    output_line = "\t".join(outline + ["1" if line_split[i] == str(a_idx) else "0" for i in include])
                                    site_store[a_idx + 1] = output_line
                                    count_store[a_idx + 1] = s_count

                    stored_counts = list(count_store.values())
                    no_vars_stored = len(stored_counts)
                    total_count = sum(stored_counts)
                    is_inverse = no_vars_stored == 2 and total_count >= no_samples * freq_inv
                    if not store_inv and is_inverse:
                        first_idx = min(count_store, key=count_store.get)
                        out_file.write(site_store[first_idx] + "\n")
                        no_stored_vars += 1
                    else:
                        for i in sorted(site_store):
                            out_file.write(site_store[i] + "\n")
                            no_stored_vars += 1

                    if no_sites % 1000 == 0:
                        print(f"\r - {no_sites} variant sites processed.")

    except Exception as e:
        print(f"VCF file did not open: {e}")
        raise typer.Exit(code=1)

    print(f"\r - {no_sites} variant sites processed.\n\nSummary:")
    print(f" - {no_sites} sites contained {no_variants} variants.")
    print(f" - {n_vars} variants were gaps.")
    print(f" - {inverse} variants were inverse alleles.")
    print(f" - {no_stored_vars}/{no_variants} variants were printed to output file [freq: {freq_low}-{freq_high}].\n\n")

if __name__ == "__main__":
    typer.run(expand)
