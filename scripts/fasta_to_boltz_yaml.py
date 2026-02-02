#!/usr/bin/env python3
"""
Description : FASTA to Boltz-2 YAML Converter
Usage:
    ~$ python fasta_to_boltz_yaml.py INPUT.fasta [-lig LIGAND.fasta] [-o OUTPUT_DIR] [--quiet]
"""

import sys
import argparse
from Bio import SeqIO
from pathlib import Path


def parse_fasta(fasta_file):
    """
    Parse a FASTA file and return a list of (id, sequence) tuples.
    Args:
        fasta_file: Path to FASTA file
    Returns:
        List of tuples: [(id, sequence), ...]
    """
    sequences = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        seq_id = record.id.split('|')[1] if '|' in record.id else record.id
        sequences.append((seq_id, str(record.seq)))

    return sequences


def write_yaml(seq, smiles, output_file):
    """
    Write Boltz-2 compatible YAML with proteins, ligand and complexes.
    Args:
        sequences: List of (id, sequence) tuples
        ligand: (id, smiles) tuple or None
        output_file: Path to output YAML file
    """
    with open(output_file, 'w') as f:
        f.write("version: 1\n")
        f.write("sequences:\n")

        # Proteins
        f.write("  - protein:\n")
        f.write(f"      id: A\n")
        f.write(f"      sequence: {seq}\n")

        # Ligand (optional)
        if smiles:
            f.write(f"  - ligand:\n")
            f.write(f"      id: L\n")
            f.write(f"      smiles: {smiles}\n")
        # Complexes (enable affinity)
            f.write("properties:\n")
            f.write(f"  - affinity:\n")
            f.write(f"      binder: L\n")
        

def fasta_to_boltz_yaml(fasta_file, ligand_file=None, output_dir=None, verbose=True):
    """
    Convert FASTA file to Boltz-2 YAML format.
    
    Args:
        fasta_file: Path to input FASTA file
        output_dir: Path to output YAML file directory (default: current)
        verbose: Print progress messages
        
    Returns:
        Path to output file
    """
    # Setup paths
    fasta_path = Path(fasta_file)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")
    
    ligand = None
    if ligand_file:
        ligs = parse_fasta(ligand_file)
        if len(ligs) != 1:
            raise ValueError("Ligand FASTA must contain exactly one entry")
        ligand = ligs[0]

    if output_dir is None:
        output_dir = Path(".")
    
    # Parse FASTA
    if verbose:
        print(f"Reading FASTA file: {fasta_path}\n")
    
    sequences = parse_fasta(fasta_path)
    
    if not sequences:
        raise ValueError("No sequences found in FASTA file")
    
    if verbose:
        print(f"Found {len(sequences)} protein sequences:")
        for i, (seq_id, seq) in enumerate(sequences, 1):
            print(f"  [{i}] {seq_id} - {len(seq)} amino acids")
        if ligand:
            print(f"Ligand: {ligand[0]}\n")
    
    # Write YAML
    for seq_id, seq in sequences:
        if ligand:
            lig_id, smiles = ligand
            yaml_file = Path(output_dir) / f"{seq_id}_{lig_id}.yaml"
        else:
            smiles = None
            yaml_file = Path(output_dir) / f"{seq_id}.yaml"
        
        if verbose:
            print(f"Writing YAML file: {yaml_file}")

        # Appel à write_yaml pour une seule séquence sous forme de liste [(seq_id, seq)]
        write_yaml(seq, smiles, yaml_file)

    if verbose:
        print(f"\nSuccessfully converted {len(sequences)} sequences!")
        print(f"  Input:  {fasta_path}")
        print(f"  Output directory: {output_dir}")


def main():
    """
    Command-line interface for FASTA to Boltz-2 YAML conversion.
    """
    parser = argparse.ArgumentParser(
        description='Convert FASTA file to Boltz-2 compatible YAML format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
--- Usage ---
  Basic conversion:
  ~$ python fasta_to_boltz_yaml.py input.fasta (-lig ligand.fasta)
  
  Specify output file:
  ~$ python fasta_to_boltz_yaml.py input.fasta (-lig ligand.fasta) -o OUTPUT_DIR/SEQ_ID.yaml
  
  Silent mode:
  ~$ python fasta_to_boltz_yaml.py input.fasta (-lig ligand.fasta) --quiet

--- Output Format ---
  
  Individual sequences found in FASTA:
    sequences:
      - protein:
          id: protein_1
          sequence: MKTAYIAK...
      - protein:
          id: protein_2
          sequence: MKTAYIAK...
          
        """
    )
    
    parser.add_argument("fasta_file", help="Protein FASTA file")
    parser.add_argument("-lig", "--ligand", help="Ligand FASTA file (SMILES)")
    parser.add_argument("-o", "--output", help="Output YAML file")
    parser.add_argument("-q", "--quiet", action="store_true", help="Silent mode")
    args = parser.parse_args()

    try:
        fasta_to_boltz_yaml(
            args.fasta_file,
            ligand_file=args.ligand,
            output_dir=args.output,
            verbose=not args.quiet
        )
        return 0
    
    except Exception as e:
        print(f"[!] Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())