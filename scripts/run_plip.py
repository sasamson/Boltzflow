#!/usr/bin/env python3
# coding: utf-8

#--------------------------------------------------------#
# Program: Run PLIP program on PDB Protein-ligand complex generated with AutoDock Tools 4
# Usage:
#   ~$ python3 run_plip.py [CPLX_FILEPATH] [OUTDIR]
#
# Author: Samantha SAMSON
# Date: 24/04/2020
#       03/01/2022 --> update docstrings
#
# Docs PLIP: https://github.com/pharmai/plip/blob/master/DOCUMENTATION.md
#--------------------------------------------------------#

# Python standard library
import sys
import os
import shutil
import subprocess


def main():
    """
    Run PLIP on a protein–ligand complex PDB file.

    Inputs:
        - cplx_pattern: identifier of the complex (for logging)
        - cplx_path: path to protein–ligand PDB file
        - outdir: output directory for PLIP results

    Outputs:
        - PLIP output files in outdir
    """

    if len(sys.argv) != 4:
        sys.stderr.write(
            "Usage: ~$ python run_plip.py [complex_id] [pdb_path] [output_dir]\n"
        )
        return 1

    cplx_id = sys.argv[1]
    cplx_path = sys.argv[2]
    outdir = sys.argv[3]

    if not os.path.isfile(cplx_path):
        sys.stderr.write(f"[Error] PDB file not found: {cplx_path}\n")
        return 1

    # Clean or create output directory
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir, exist_ok=True)

    print(f"\n>>> PLIP interaction profiling: {cplx_id}\n")
    print("Command: plip")
    print("Options: -f [PDB_FILE] -xty -o [OUTPUT_DIR]")
    print(f"> FILE: {cplx_path}")
    print(f"> OUTDIR: {outdir}\n")

    cmd = [
        "plip",
        "-f", cplx_path,
        "-xty",
        "-o", outdir
    ]

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"PLIP failed with error:\n{e}\n")
        return 1

    print("PLIP analysis completed successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())