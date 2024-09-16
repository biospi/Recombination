import pandas as pd
import numpy as np
from pathlib import Path
import plotly.graph_objs as go
import plotly.io as pio
import plotly.express as px
import ast
from utils import run_cmd
import matplotlib.pyplot as plt


def convert_to_dict(count_str):
    count_str = count_str.replace("  ", ", ")
    return ast.literal_eval(count_str)


def histogram_fasta(input_file):
    df = pd.read_csv(input_file)
    # Data for plotting
    x = df['record_id']
    df['count'] = df['count'].apply(convert_to_dict)
    counts = df['count'].apply(lambda x: dict(x).get('N', 0))


    plt.figure(figsize=(10, 6))
    plt.hist(counts, bins=1000)
    plt.title('Histogram of N in fasta')
    plt.xlabel('Number of Ns')
    plt.ylabel('Frequency')
    # plt.show()
    output_file = input_file.parent / f"hisrogram_fasta.png"
    print(output_file)
    plt.savefig(output_file) 
    plt.close()
    
    # # Calculate the frequency of each count
    # count_frequency = counts.value_counts().sort_index()

    # # Create a bar chart
    # fig = go.Figure(data=[go.Bar(x=count_frequency.index, y=count_frequency.values)])

    # # Update layout
    # fig.update_layout(
    #     title='Frequency of Missing Values (N) per Count',
    #     xaxis_title='Count of Missing Values (N)',
    #     yaxis_title='Frequency',
    #     xaxis_tickangle=-45
    # )
    # fig.update_traces(dict(marker_line_width=0))

    # # Save the figure as an HTML file
    # output_html = input_file.parent / 'missing_values_plot.html'
    # print(output_html)
    # pio.write_html(fig, file=output_html)


def histogram_vcf(input_file, title="title"):
    print(input_file)
    missing_snp_counts = []

    with open(input_file, 'r') as file:
        cpt = 0
        tot = 0
        iso_list = []
        for i, line in enumerate(file):
            if line.startswith('##'):
                continue

            if line.startswith('#CHROM'):
                columns = line.strip().split('\t')
                sample_names = columns[9:]  # Sample names start from the 10th column onward
                missing_snp_counts = [0] * len(sample_names) 
            #print('\t'.join([x[0:4] for x in line.split('\t')][0:27]))  

            tot += 1
            columns = line.strip().split('\t')
            genotypes = columns[9:]  # Genotype data starts from the 10th column onward
            
            # Check each genotype for missing data (represented as './.')
            if "*" in line:
                for i, genotype in enumerate(genotypes):
                    if genotype == '1':  # Check for missing SNPs
                        missing_snp_counts[i] += 1
                        #print(missing_snp_counts)
    
    plt.figure(figsize=(10, 6))
    plt.hist(missing_snp_counts, bins=100)
    plt.title(f'Histogram of Missing SNPs\n{title}')
    plt.xlabel('Number of Missing SNPs')
    plt.ylabel('Frequency')
    output_file = input_file.parent / f"hisrogram_{title}.png"
    print(output_file)
    plt.savefig(output_file) 
    plt.close()


    #     print(f"Total after filtration = {cpt}/{tot}")
    #     df = pd.DataFrame(iso_list, columns=["isolate"])

    #     count_frequency = df["isolate"][1:].value_counts().sort_index()

    #     # Create a bar chart
    #     fig = go.Figure(data=[go.Bar(x=count_frequency.index, y=count_frequency.values)])

    #     # Update layout
    #     fig.update_layout(
    #         title='Frequency of Missing Values (*) in VCF',
    #         xaxis_title='Isolate with (*)',
    #         yaxis_title='Frequency',
    #         xaxis_tickangle=-45
    #     )
    #     fig.update_traces(dict(marker_line_width=0))

    #     # Save the figure as an HTML file
    #     output_html = input_file.parent / 'missing_values_plot_vcf.html'
    #     print(output_html)
    #     pio.write_html(fig, file=output_html)


        # fig = px.histogram(df, x="isolate", barmode="overlay", title="Histogram of N_MISS and F_MISS")
        # fig.update_layout(
        #     title='missingness (*) in vcf'.title(),
        #     xaxis_tickangle=-45
        # )
        # fig.update_traces(dict(marker_line_width=0))
        # output_html = input_file.parent / 'vcf_hist.html'
        # print(output_html)
        # pio.write_html(fig, file=output_html)


    # if not Path("output.imiss").exists():
    # run_cmd(f"vcftools --vcf {input_file.as_posix()} --out {input_file.parent.as_posix()} --missing-indv", input_file.parent,"vcf", -1)
    # run_cmd(f"vcftools --vcf {input_file.as_posix()} --out {input_file.parent.as_posix()} --missing-site", input_file.parent,"vcf", -1)

    # df = pd.read_csv("output.lmiss", sep='\t')
    
    # x = df['POS']
    # y = df['N_MISS']
    # fig = go.Figure(data=[go.Bar(x=x, y=y)])
    # fig.update_layout(
    #     title='missingness on a per-individual basis'.title(),
    #     xaxis_title='Record ID',
    #     yaxis_title='N_MISS',
    #     xaxis_tickangle=-45
    # )
    # fig.update_traces(dict(marker_line_width=0))
    # output_html = input_file.parent / 'imiss.html'
    # print(output_html)
    # pio.write_html(fig, file=output_html)
    

def count_missing_snps(vcf_file):
    # Initialize variables
    sample_names = []
    missing_snp_counts = []

    # Open and read the VCF file line by line
    with open(vcf_file, 'r') as file:
        for line in file:
            # Skip header lines (they start with ##)
            if line.startswith('##'):
                continue

            # Process the header line to get sample names
            if line.startswith('#CHROM'):
                columns = line.strip().split('\t')
                sample_names = columns[9:]  # Sample names start from the 10th column onward
                missing_snp_counts = [0] * len(sample_names)  # Initialize a list for missing counts
                continue
            
            # For all other lines, process the SNP data
            columns = line.strip().split('\t')
            genotypes = columns[9:]  # Genotype data starts from the 10th column onward
            
            # Check each genotype for missing data (represented as './.')
            for i, genotype in enumerate(genotypes):
                if genotype.startswith('./.') or genotype == '.' or genotype == 'N':  # Check for missing SNPs
                    missing_snp_counts[i] += 1

    return sample_names, missing_snp_counts


def plot_histogram(missing_snp_counts):
    a = np.array(missing_snp_counts)
    a = a[a>0]
    print(a)
    # Create a histogram of missing SNPs
    plt.figure(figsize=(10, 6))
    plt.hist(a, bins=1)
    plt.title('Histogram of Missing SNPs')
    plt.xlabel('Number of Missing SNPs')
    plt.ylabel('Frequency')
    plt.show()

if __name__ == "__main__":
    histogram_fasta(Path('/home/axel/python-projects/Recombination/output/log_4.csv'))
    #histogram_vcf(Path('output/wgs-mapping_with_ref.fasta.expand.vcf'))