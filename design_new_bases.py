import os
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem import AllChem
from ase.io import write
from ase import Atoms

# Ensure RDKit can generate 3D coordinates
rdDepictor.SetPreferCoordGen(True)


# Define function to generate non-natural nucleotide
def generate_non_natural_nucleotide():
    # Create a base with non-natural modification
    base_smiles = "c1ccccc1"  # Benzene ring as a placeholder for a non-natural base
    base = Chem.MolFromSmiles(base_smiles)

    # Add a non-natural sugar
    sugar_smiles = "C1COC(O1)(C2CC2)"  # Placeholder for modified sugar
    sugar = Chem.MolFromSmiles(sugar_smiles)

    # Add a non-standard phosphate group
    phosphate_smiles = "OP(=O)(O)OC3CC3"  # Phosphorothioate group as placeholder
    phosphate = Chem.MolFromSmiles(phosphate_smiles)

    # Combine the base, sugar, and phosphate into one molecule
    nucleotide = Chem.CombineMols(base, sugar)
    nucleotide = Chem.CombineMols(nucleotide, phosphate)

    # Add artificial functional group (e.g., fluorine atom)
    functional_group_smiles = "F"
    functional_group = Chem.MolFromSmiles(functional_group_smiles)
    nucleotide = Chem.CombineMols(nucleotide, functional_group)

    # Generate 3D coordinates
    nucleotide = Chem.AddHs(nucleotide)
    AllChem.EmbedMolecule(nucleotide, randomSeed=42)
    AllChem.UFFOptimizeMolecule(nucleotide)

    return nucleotide


# Generate and save non-natural nucleotides
def save_nucleotides(n, prefix="nucleotide"):
    if not os.path.exists("output"):
        os.makedirs("output")

    for i in range(n):
        nucleotide = generate_non_natural_nucleotide()
        atoms = Atoms(Chem.MolToXYZBlock(nucleotide))
        file_name = f"output/{prefix}_{i + 1}.gjf"
        write(file_name, atoms)


# Generate dimers, trimers, or tetramers of nucleotides
def generate_oligomers(nucleotide, oligomer_size):
    mol = nucleotide
    for _ in range(oligomer_size - 1):
        mol = Chem.CombineMols(mol, nucleotide)
    return mol


def save_oligomers(n, oligomer_size, prefix="oligomer"):
    if not os.path.exists("output"):
        os.makedirs("output")

    for i in range(n):
        nucleotide = generate_non_natural_nucleotide()
        oligomer = generate_oligomers(nucleotide, oligomer_size)
        atoms = Atoms(Chem.MolToXYZBlock(oligomer))
        file_name = f"output/{prefix}_{oligomer_size}_{i + 1}.gjf"
        write(file_name, atoms)


# Example usage
save_nucleotides(10)  # Generate 10 non-natural nucleotides
save_oligomers(5, 2)  # Generate 5 dimers
save_oligomers(5, 3)  # Generate 5 trimers
save_oligomers(5, 4)  # Generate 5 tetramers
