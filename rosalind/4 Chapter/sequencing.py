"""
#! Cyclopeptide Sequencing Problem

Given an ideal experimental spectrum, find a cyclic peptide whose theoretical spectrum matches the experimental spectrum.

Given: A collection of (possibly repeated) integers Spectrum corresponding to an ideal experimental spectrum.
Return: Every amino acid string Peptide such that Cyclospectrum(Peptide) = Spectrum (if such a string exists).
"""
from typing import List

INTEGER_MASS = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
                'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}


def expand_peptides(peptides: List[str], expander=INTEGER_MASS.keys()):
    """
    Extends each peptide in given list by 1 amino acid
    @param: peptides: List[str] -- given list of peptides
    @param: expander: Dict_keys -- elements by which we expand
    """
    res = []
    for pep in peptides:
        for acid in expander:
            res.append(pep + acid)
    return res


def calculate_spectrum(peptide: str):
    """
    Calculates theoretical spectrum for provided peptide
    @param: peptide: str -- provided peptide
    """
    n = len(peptide)
    cyclic_peptide = peptide+peptide
    spectrum = [0] * (n * (n-1))
    cyclic_subpeptides = [None] * (n * (n-1))
    ind = 0
    for i in range(n):
        for j in range(n - 1):
            subpeptide = cyclic_peptide[i: i+j+1]
            spectrum[ind] = sum([INTEGER_MASS[acid]
                                 for acid in subpeptide])
            cyclic_subpeptides[ind] = subpeptide
            ind += 1
    spectrum.append(sum([INTEGER_MASS[acid] for acid in peptide]))
    spectrum.append(0)
    cyclic_subpeptides.extend([' ', peptide])
    return spectrum, cyclic_subpeptides


def linear_spectrum(peptide: str, verbose=False):
    """
    Calculates linear spectrum for given peptide
    @param: peptide: str

    LinearSpectrum(Peptide, Alphabet, AminoAcidMass)
        PrefixMass(0) ← 0
        for i ← 1 to |Peptide|
            for every symbol s in Alphabet
                if s = i-th amino acid in Peptide
                    PrefixMass(i) ← PrefixMass(i − 1) + AminoAcidMass[s]
        LinearSpectrum ← a list consisting of the single integer 0
        for i ← 0 to |Peptide| − 1
            for j ← i + 1 to |Peptide|
                add PrefixMass(j) − PrefixMass(i) to LinearSpectrum
        return sorted list LinearSpectrum
    """
    n = len(peptide)
    prefix_mass = [0] * (n+1)
    for i in range(1, n+1):
        acid = peptide[i-1]
        prefix_mass[i] = prefix_mass[i-1] + INTEGER_MASS[acid]
    if verbose:
        print(prefix_mass)
    linear_spec = [0] * int((n * (n+1)/2) + 1)
    ind = 0
    for i in range(0, n):
        for j in range(i+1, n+1):
            linear_spec[ind] = prefix_mass[j] - prefix_mass[i]
            ind += 1
    return sorted(linear_spec)


def cyclopeptide_sequencing(spectrum: List[int]):
    """
    Finds a cyclic peptide whose theoretical spectrum matches the experimental spectrum
    @param: spectrum: List[int] -- given experimental spectrum

     CyclopeptideSequencing(Spectrum)
        CandidatePeptides ← a set containing only the empty peptide
        FinalPeptides ← empty list of strings
        while CandidatePeptides is nonempty
            CandidatePeptides ← Expand(CandidatePeptides)
            for each peptide Peptide in CandidatePeptides
                if Mass(Peptide) = ParentMass(Spectrum)
                    if Cyclospectrum(Peptide) = Spectrum and Peptide is not in FinalPeptides
                        append Peptide to FinalPeptides
                    remove Peptide from CandidatePeptides
                else if Peptide is not consistent with Spectrum
                    remove Peptide from CandidatePeptides
        return FinalPeptides
    """
    candidates = [""]
    finals = []
    expander = [key for key in INTEGER_MASS.keys() if key not in ["L", "Q"]]
    while candidates != []:
        candidates = expand_peptides(candidates, expander)
        candidates_to_remove = []
        for pep in candidates:
            spec, _ = calculate_spectrum(pep)
            linear_spec = linear_spectrum(pep)
            if linear_spec[-1] == spectrum[-1]:
                if sorted(spec) == sorted(spectrum):
                    finals.append(pep)
                candidates_to_remove.append(pep)
            elif not all(lin in spectrum for lin in linear_spec):
                candidates_to_remove.append(pep)
        candidates = [
            pep for pep in candidates if pep not in candidates_to_remove]
        print(candidates)
    print(finals)
    return finals


def main():
    with open("/home/snopoff/Downloads/rosalind_ba4e.txt", "r") as f:
        spectrum = [int(number) for number in f.readline().strip().split(' ')]
    peptides = cyclopeptide_sequencing(spectrum)
    masses = [[INTEGER_MASS[acid] for acid in pep] for pep in peptides]
    formatted_masses = ["-".join(list(map(str, m))) for m in masses]
    res = " ".join(formatted_masses)
    with open("res.txt", "w") as f:
        f.write(" ".join(sorted(res.split(" "))))


if __name__ == "__main__":
    main()
