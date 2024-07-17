from pathlib import Path
import typer
from typing import List, Optional
import warnings

from utils import count_isolates_in_tar_gz, filter_isolates_in_tar_gz, run_cmd
warnings.filterwarnings("ignore")
import pandas as pd


def main(
    dataset_filepath: Path = typer.Option(
        ..., exists=True, file_okay=True, dir_okay=False, resolve_path=True
    ),
    output_dir: str = "output",
    exlude: List[str] = ["SRR5193868", "SRR5193995", "SRR7850489", "SRR5193283", "SRR3049562", "SRR7495458"], 
):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"dataset_filepath={dataset_filepath}")
    print("loading...")
    iso_count = count_isolates_in_tar_gz(dataset_filepath.as_posix())
    print(f"Found {iso_count} isolates in {dataset_filepath}.")

    cleaned_tar_gz_filepath = dataset_filepath.parent / "cleaned-wgs-mapping.tar.gz"
    if cleaned_tar_gz_filepath.exists():
        print(f"Deleting {cleaned_tar_gz_filepath}...")
        cleaned_tar_gz_filepath.unlink()

    filter_isolates_in_tar_gz(dataset_filepath.as_posix(), cleaned_tar_gz_filepath.as_posix(), exlude)
    
    iso_count = count_isolates_in_tar_gz(cleaned_tar_gz_filepath.as_posix())
    print(f"Found {iso_count} isolates in {cleaned_tar_gz_filepath}.")

    filepath_vcf = Path(f"{cleaned_tar_gz_filepath.name}.vcf")
    print(filepath_vcf)
    print("building snp file...")
    run_cmd(
        f"snp-sites -o {filepath_vcf.as_posix()} -v {cleaned_tar_gz_filepath.as_posix()}",
        output_dir,
        "vcf",
        iso_count,
    )


if __name__ == "__main__":
    #dataset_filepath = Path("wgs-mapping.tar.gz")
    #main(dataset_filepath)
    typer.run(main)
