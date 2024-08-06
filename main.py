from pathlib import Path
import typer
from typing import List, Optional
import warnings

from cluster import cluster
from utils import count_isolates_in_tar_gz, filter_isolates_in_tar_gz, get_header, run_cmd
warnings.filterwarnings("ignore")
import pandas as pd


def main(
    dataset_filepath: Path = typer.Option(
        ..., exists=True, file_okay=True, dir_okay=False, resolve_path=True
    ),
    output_dir: str = "output",
    # exlude: List[str] = ["SRR5193868", "SRR5193995", "SRR7850489", "SRR5193283", "SRR3049562", "SRR7495458", "SRR8437386", "SRR7828233", "SRR7533537", "SRR5216167", "SRR1965543"], 
):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"dataset_filepath={dataset_filepath}")
    print("loading...")
    iso_count = count_isolates_in_tar_gz(dataset_filepath.as_posix())
    print(f"Found {iso_count} isolates in {dataset_filepath}.")

    get_header(dataset_filepath.as_posix())
    isolates, cluster_list = cluster()

    cleaned_filepath = dataset_filepath.parent / "cleaned-wgs-mapping.fasta"
    if cleaned_filepath.exists():
        print(f"Deleting {cleaned_filepath}...")
        cleaned_filepath.unlink()

    iso_count = filter_isolates_in_tar_gz(dataset_filepath.as_posix(), cleaned_filepath.as_posix(), isolates)
    
    # # iso_count = count_isolates_in_tar_gz(cleaned_tar_gz_filepath.as_posix())
    # # print(f"Found {iso_count} isolates in {cleaned_tar_gz_filepath}.")


    filepath_vcf = Path(f"{dataset_filepath.name}.vcf")
    print(filepath_vcf)
    print("building snp file...")
    run_cmd(
        f"snp-sites -o {filepath_vcf.as_posix()} -v {cleaned_filepath.as_posix()}",
        output_dir,
        "vcf",
        iso_count,
    )


if __name__ == "__main__":
    dataset_filepath = Path("wgs-mapping.tar.gz")
    #dataset_filepath = Path("mapping-NC_011294.snp_sites.fasta")
    main(dataset_filepath)
    #typer.run(main)
