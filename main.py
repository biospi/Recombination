from pathlib import Path
import typer
from typing import List, Optional
import warnings

from cluster import cluster
from utils import count_isolates_in_tar_gz, filter_isolates_in_tar_gz, get_header, run_cmd
from vcf2cp import vcf2cp
warnings.filterwarnings("ignore")
import pandas as pd


def main(
    dataset_filepath: Path = typer.Option(
        ..., exists=True, file_okay=True, dir_okay=False, resolve_path=True
    ),
    output_dir: str = "output"
):
    filepath_vcf = Path(f"{dataset_filepath.name}.vcf")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if not filepath_vcf.exists():
        print(f"dataset_filepath={dataset_filepath}")
        print("loading...")
        iso_count = count_isolates_in_tar_gz(dataset_filepath)
        print(f"Found {iso_count} isolates in {dataset_filepath}.")

        get_header(output_dir, dataset_filepath.as_posix())
        isolates, cluster_list = cluster(output_dir)

        cleaned_filepath = dataset_filepath.parent / "cleaned-wgs-mapping.fasta"
        if cleaned_filepath.exists():
            print(f"Deleting {cleaned_filepath}...")
            cleaned_filepath.unlink()

        iso_count = filter_isolates_in_tar_gz(output_dir, dataset_filepath.as_posix(), cleaned_filepath.as_posix(), isolates)
        
        # # iso_count = count_isolates_in_tar_gz(cleaned_tar_gz_filepath.as_posix())
        # # print(f"Found {iso_count} isolates in {cleaned_tar_gz_filepath}.")
        print(filepath_vcf)
        print("building snp file...")
        run_cmd(
            f"snp-sites -o {filepath_vcf.as_posix()} -v {cleaned_filepath.as_posix()}",
            output_dir,
            "vcf",
            iso_count,
        )

    print("vcf2cp...")
    vcf2cp(input_vcf=filepath_vcf.as_posix(), ploidy=1, output_prefix=(output_dir / "vcf2cp_").as_posix())


if __name__ == "__main__":
    # dataset_filepath = Path("wgs-mapping.tar.gz")
    # #dataset_filepath = Path("mapping-NC_011294.snp_sites.fasta")
    # main(dataset_filepath)
    typer.run(main)
