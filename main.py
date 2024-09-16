from pathlib import Path
import typer
from typing import List, Optional
import warnings

from cluster import cluster
import expand_vcf
from makeuniformrecfile import makeuniformrecfile
from utils import count_isolates_in_tar_gz, filter_isolates_in_tar_gz, filter_vcf, get_header, run_cmd
from vcf2cp import vcf2cp
warnings.filterwarnings("ignore")
import pandas as pd
import time
import datetime


def main(
    dataset_filepath: Path = typer.Option(
        ..., exists=True, file_okay=True, dir_okay=False, resolve_path=True
    ),
    output_dir: str = "output2"
):

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filepath_vcf = output_dir / f"{dataset_filepath.name}.vcf"

    if not filepath_vcf.exists():
        print(f"dataset_filepath={dataset_filepath}")
        print("loading...")
        iso_count = -1
        # iso_count = count_isolates_in_tar_gz(dataset_filepath)
        # print(f"Found {iso_count} isolates in {dataset_filepath}.")

        # get_header(output_dir, dataset_filepath.as_posix())
        # isolates, cluster_list = cluster(output_dir)

        # cleaned_filepath = dataset_filepath.parent / "wgs-mapping.fasta"
        # if cleaned_filepath.exists():
        #     print(f"Deleting {cleaned_filepath}...")
        #     cleaned_filepath.unlink()

        # iso_count = filter_isolates_in_tar_gz(output_dir, dataset_filepath.as_posix(), cleaned_filepath.as_posix(), isolates)
        
        # # iso_count = count_isolates_in_tar_gz(cleaned_tar_gz_filepath.as_posix())
        # # print(f"Found {iso_count} isolates in {cleaned_tar_gz_filepath}.")
        print(filepath_vcf)
        print("building snp file...")
        start = time.time()
        run_cmd(
            f"snp-sites -o {filepath_vcf.as_posix()} -v {dataset_filepath.as_posix()}",
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


    filepath_vcf_filtered = output_dir / f"{dataset_filepath.name}.filtered.vcf"    
    filter_vcf(filepath_vcf, filepath_vcf_filtered)

    filepath_vcf_expand = output_dir / f"{dataset_filepath.name}.expand.vcf"
    if not filepath_vcf_expand.exists():
        expand_vcf.expand(filepath_vcf_filtered.as_posix(),
                        filepath_vcf_expand.as_posix(),
                        freq_low=0.0,
                        freq_high=1.0,
                        n_sites=False,
                        store_inv=False,
                        sample_list=None,
                        freq_inv=0.98,
                        data_type="auto")

    print("vcf2cp...")
    start = time.time()
    vcf2cp(input_vcf=filepath_vcf_expand.as_posix(), 
           output_prefix=(output_dir / "vcf2cp_").as_posix(),
           jitter=True,
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

    print("finestructure...")
    start = time.time()
    ids = (output_dir / "vcf2cp_.ids").as_posix()
    recomb = (output_dir / "vcf2cp_.recomb").as_posix()
    phase = (output_dir / "vcf2cp_.phase").as_posix()
    out = output_dir.as_posix()
    #chmod 777 ./fs_4.1.1/fs_linux_glibc2.3
    run_cmd(
        f"./fs_linux_glibc2.3 cp -j -b -t {ids} -a 1 1 -r {recomb} -n 0.001 -M 0.001 -g {phase} -o {out} -v",
        output_dir,
        "fs"
    )
    end = time.time()
    print(
        "time (H:M:S.mm)= "
        + str(datetime.timedelta(seconds=end - start))
    )


if __name__ == "__main__":
    dataset_filepath = Path("output/wgs-mapping_with_ref.fasta")
    main(dataset_filepath, output_dir="output")
    #typer.run(main)
