import sys
from pathlib import Path

def parse_vcf(vcf_file):
    """Parse a VCF file and return its header and data."""
    with open(vcf_file, 'r') as file:
        header = []
        data = []
        for line in file:
            if line.startswith("#"):  # Header lines
                header.append(line.strip())
            else:  # Data lines
                data.append(line.strip().split("\t"))
    return header, data

def calculate_missingness(site):
    """Calculate the percentage of missing data (represented by '*') in a site."""
    genotypes = site[9:]  # All genotype fields (after FORMAT)
    total_genotypes = len(genotypes)
    missing_count = sum(1 for gt in genotypes if '*' in gt)  # Count '*' alleles
    missingness = (missing_count / total_genotypes) * 100
    return missingness

def filter_sites_by_missingness(data, threshold=10.0):
    """Filter sites by missingness threshold."""
    filtered_data = []
    for site in data:
        missingness = calculate_missingness(site)
        if missingness <= threshold:
            filtered_data.append(site)
    return filtered_data

def write_filtered_vcf(output_file, header, filtered_data):
    """Write the filtered VCF data to a new file."""
    with open(output_file, 'w') as file:
        # Write header
        for line in header:
            file.write(line + "\n")
        # Write filtered data
        for site in filtered_data:
            file.write("\t".join(site) + "\n")

def main(input_vcf, output_vcf, threshold=10.0):
    # Parse the VCF file
    header, data = parse_vcf(input_vcf)
    
    # Filter sites based on missingness
    filtered_data = filter_sites_by_missingness(data, threshold)
    
    # Write the filtered VCF to a new file
    write_filtered_vcf(output_vcf, header, filtered_data)
    print(f"Filtered VCF written to {output_vcf}")

if __name__ == "__main__":
    input_vcf = Path('output/wgs-mapping_with_ref.fasta.vcf')
    output_vcf = Path('output/wgs-mapping_with_ref.filtered.fasta.vcf')
    main(input_vcf, output_vcf)