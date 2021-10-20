"""
#! Counting Peptides with Given Mass Problem
Denote the total mass of an amino acid string Peptide as Mass(Peptide). 
In mass spectrometry experiments, whereas the peptide that generated a spectrum is unknown, 
the peptide’s mass is typically known and is denoted ParentMass(Spectrum). 
Of course, given an ideal experimental spectrum, Mass(Peptide) is given by the largest mass in the spectrum.
A brute force approach to reconstructing a peptide from its theoretical spectrum would generate all possible peptides 
whose mass is equal to ParentMass(Spectrum) and then check which of these peptides has theoretical spectra matching Spectrum. 
However, we should be concerned about the running time of such an approach: 
how many peptides are there having mass equal to ParentMass(Spectrum)?

Compute the number of peptides of given total mass.

Given: An integer m.
Return: The number of linear peptides having integer mass m.
"""
from collections import defaultdict
from typing import DefaultDict

INTEGER_MASS = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
                'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}


def bf_cyclopeptide_sequencing(mass: int, pos_mass_dict: DefaultDict[int, int]):
    """
    Bruteforcely computes the number of peptides of given mass
    @param: mass: int -- given mass
    @param: pos_mass_dict: DefaultDict[int] -- provided dict of possible masses and the amount of linear peptides

    BFCyclopeptideSequencing(Spectrum)
        for every peptide with Mass(Peptide) equal to ParentMass(Spectrum)
            if Spectrum = Cyclospectrum(Peptide)
                output Peptide
    """
    if mass == 0:
        return 1, pos_mass_dict
    if mass < 57:
        return 0, pos_mass_dict
    if mass in pos_mass_dict:
        return pos_mass_dict[mass], pos_mass_dict
    count = 0
    for weight in (set(INTEGER_MASS.values())):
        amount, pos_mass_dict = bf_cyclopeptide_sequencing(
            mass - weight, pos_mass_dict)
        count += amount
    pos_mass_dict[mass] = count
    return count, pos_mass_dict


def main():
    with open("/home/snopoff/Downloads/rosalind_ba4d.txt", "r") as f:
        mass = int(f.readline().strip())
    pos_mass_dict = defaultdict(int)
    amount, _ = bf_cyclopeptide_sequencing(
        mass, pos_mass_dict)
    with open("res.txt", "w") as f:
        f.write(str(amount))


if __name__ == "__main__":
    main()
