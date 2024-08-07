import typer
import sys
from pathlib import Path
from tqdm import tqdm


def trim(string):
    """Remove whitespace from the beginning and end of the argument."""
    return string.strip()

def help():
    """Print the help message and exit."""
    print("Usage: python makeuniformrecfile.py [OPTIONS] PHASEFILE OUTPUTFILE")
    print("     PHASEFILE: A valid chromopainter inputfile ending in .phase (in ChromoPainter v1 or v2 format)")
    print("     OUTPUTFILE: A recombination file usable with PHASEFILE in ChromoPainter, nominally in Morgans/base.")
    print("Options:")
    print("  -c, --centimorgansperbase <FLOAT>: Recombination rate per base, default is 1/10000000")
    print("The recombination rate is scaled to be approximately that in humans (1 centi-Morgan/Mb, or 100Mb for one Morgan).")
    print("Because of this, it will NOT be usable directly and should only ever be used in conjunction with EM parameter estimation, which corrects for the global amount of recombination.")
    print("If you are working on non-humans or simulated data, you may experience problems with EM estimation.")
    print("The parameter may get stuck at a local mode where there is effectively infinite (or no) recombination.")
    print("In this case, you should specify the initial conditions of ChromoPainter to have a much smaller or larger Ne (-n) value.")
    raise typer.Exit()


def makeuniformrecfile(phasefile: str, outputfile: str, centimorgansperbase: float = 1.0 / 10000000):
    print(f"phasefile={phasefile}")
    print(f"outputfile={outputfile}")

    if Path(outputfile).exists():
        print(f"output {outputfile} already exist.")
        return
    
    snplocs = ""
    snps = []
    # Read in the SNPS
    try:
        with open(phasefile, 'r') as phase_file:
            for _ in range(4):
                snplocs = trim(phase_file.readline())
                if snplocs.startswith("P"):
                    break
    except Exception as e:
        print(f"Error reading phasefile: {e}", err=True)

    snps = snplocs.split()
    numsnps = len(snps) - 1

    print(f"Found {numsnps} SNPS")
    finalSNP = snps[numsnps]

    # Set up the outputfile
    try:
        with open(outputfile, 'w') as output_file:
            output_file.write("start.pos recom.rate.perbp\n")

            # Loop over all SNPs to calculate their recombination position
            for i in tqdm(range(numsnps - 1)):
                snploc = int(snps[i + 1])
                overallrecrate = centimorgansperbase
                if i < numsnps - 1 and int(snps[i + 2]) < snploc:
                    overallrecrate = -9
                output_file.write(f"{snploc} {overallrecrate:.20f}\n")

            output_file.write(f"{finalSNP} 0\n")
    except Exception as e:
        print(f"Error writing to outputfile: {e}", err=True)

    print("Note: No absolute recombination rate.")
    print("You must perform EM estimation of Ne (e.g. -in -i 10) using chromopainter to use this recombination file")

if __name__ == "__main__":
    typer.run(makeuniformrecfile)
