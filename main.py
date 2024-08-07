from pathlib import Path
import typer
from typing import List, Optional
import warnings

from cluster import cluster
from makeuniformrecfile import makeuniformrecfile
from utils import count_isolates_in_tar_gz, filter_isolates_in_tar_gz, get_header, run_cmd
from vcf2cp import vcf2cp
warnings.filterwarnings("ignore")
import pandas as pd
import time
import datetime


def main(
    dataset_filepath: Path = typer.Option(
        ..., exists=True, file_okay=True, dir_okay=False, resolve_path=True
    ),
    output_dir: str = "output"
):

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filepath_vcf = output_dir / f"{dataset_filepath.name}.vcf"

    if not filepath_vcf.exists():
        print(f"dataset_filepath={dataset_filepath}")
        print("loading...")
        iso_count = count_isolates_in_tar_gz(dataset_filepath)
        print(f"Found {iso_count} isolates in {dataset_filepath}.")

        get_header(output_dir, dataset_filepath.as_posix())
        isolates, cluster_list = cluster(output_dir)

        cleaned_filepath = dataset_filepath.parent / "wgs-mapping.fasta"
        if cleaned_filepath.exists():
            print(f"Deleting {cleaned_filepath}...")
            cleaned_filepath.unlink()

        iso_count = filter_isolates_in_tar_gz(output_dir, dataset_filepath.as_posix(), cleaned_filepath.as_posix(), isolates)
        
        # # iso_count = count_isolates_in_tar_gz(cleaned_tar_gz_filepath.as_posix())
        # # print(f"Found {iso_count} isolates in {cleaned_tar_gz_filepath}.")
        print(filepath_vcf)
        print("building snp file...")
        start = time.time()
        run_cmd(
            f"snp-sites -o {filepath_vcf.as_posix()} -v {cleaned_filepath.as_posix()}",
            output_dir,
            "vcf",
            iso_count,
        )
        end = time.time()
        print(
            "time (H:M:S.mm)= "
            + str(datetime.timedelta(seconds=end - start))
        )
    else:
        print(f"filepath_vcf={filepath_vcf}")

    print("vcf2cp...")
    start = time.time()
    vcf2cp(input_vcf=filepath_vcf.as_posix(), 
           output_prefix=(output_dir / "vcf2cp_").as_posix(),
           jitter=False,
           numindsmaxsize=500,
           chromosomegap=-1,
           ploidy=1)
    end = time.time()
    print(
        "time (H:M:S.mm)= "
        + str(datetime.timedelta(seconds=end - start))
    )

    print("makeuniformrecfile...")
    start = time.time()
    makeuniformrecfile(phasefile=(output_dir / "vcf2cp_.phase").as_posix(),
                        outputfile=(output_dir / "vcf2cp_.recomb").as_posix()
                        )
    end = time.time()
    print(
        "time (H:M:S.mm)= "
        + str(datetime.timedelta(seconds=end - start))
    )


if __name__ == "__main__":
    #dataset_filepath = Path("wgs-mapping.tar.gz")
    #main(dataset_filepath)
    typer.run(main)
