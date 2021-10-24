"""
#! Implement LeaderboardCyclopeptideSequencing

Given: An integer N and a collection of integers Spectrum.
Return: LeaderPeptide after running LeaderboardCyclopeptideSequencing(Spectrum, N).
"""
from typing import List

INTEGER_MASS = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
                'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}


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


def cyclic_spectrum(peptide: str, verbose=False):
    """
    Calculates cyclic spectrum for provided peptide
    @param: peptide: str -- provided peptide

    CyclicSpectrum(Peptide, Alphabet, AminoAcidMass)
        PrefixMass(0) ← 0
        for i ← 1 to |Peptide|
            for every symbol s in Alphabet
                if s = i-th amino acid in Peptide
                    PrefixMass(i) ← PrefixMass(i − 1) + AminoAcidMass﻿[s]
        peptideMass ← PrefixMass(|Peptide|)
        CyclicSpectrum ← a list consisting of the single integer 0
        for i ← 0 to |Peptide| − 1
            for j ← i + 1 to |Peptide|
                add PrefixMass(j) − PrefixMass(i) to CyclicSpectrum
                if i > 0 and j < |Peptide|
                    add peptideMass - (PrefixMass(j) − PrefixMass(i)) to CyclicSpectrum
        return sorted list CyclicSpectrum
    """
    n = len(peptide)
    prefix_mass = [0] * (n+1)
    for i in range(1, n+1):
        acid = peptide[i-1]
        prefix_mass[i] = prefix_mass[i-1] + INTEGER_MASS[acid]
    if verbose:
        print(prefix_mass)
    cyclic_spec = [0]
    for i in range(0, n):
        for j in range(i+1, n+1):
            cyclic_spec.append(prefix_mass[j] - prefix_mass[i])
            if i > 0 and j < n:
                cyclic_spec.append(
                    prefix_mass[-1] - (prefix_mass[j] - prefix_mass[i]))
    return sorted(cyclic_spec)


def compute_score(peptide: str, spectrum: List[int]):
    """
    Computes score for given peptide and spectrum
    @param: peptide: str -- given peptide
    @param: spectrum: List[str] -- given spectrum
    """
    spec = cyclic_spectrum(peptide)
    sc = 0
    for s in spectrum:
        if s in spec:
            sc += s in spec
            spec.remove(s)
    return sc


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


def trim(leaderboard: List[str], spectrum: List[int], N: int):
    """
    Cut the leaderboard up to first best N score values
    @param: leaderboard: List[str] -- given list of peptide
    @param: spectrum: List[int] -- given spectrum
    @param: N: int -- integer value which we consider to cut to
    """
    scores = sorted([compute_score(peptide, spectrum)
                    for peptide in leaderboard], reverse=True)
    if N > len(scores):
        return leaderboard
    lowest_score = scores[N]
    result = []
    for peptide in leaderboard:
        if compute_score(peptide, spectrum) >= lowest_score:
            result.append(peptide)
    return result


def leaderboard_cyclopeptide_sequencing(spectrum: List[int], N: int):
    """
    Finds a cyclic peptide whose theoretical spectrum matches the experimental spectrum 
    with so called Leaderboard 
    @param: spectrum: List[int] -- given spectrum 
    @param: N: int -- number of best N scores we consider

    LeaderboardCyclopeptideSequencing(Spectrum, N)
        Leaderboard ← set containing only the empty peptide
        LeaderPeptide ← empty peptide
        while Leaderboard is non-empty
            Leaderboard ← Expand(Leaderboard)
            for each Peptide in Leaderboard
                if Mass(Peptide) = ParentMass(Spectrum)
                    if Score(Peptide, Spectrum) > Score(LeaderPeptide, Spectrum)
                        LeaderPeptide ← Peptide
                else if Mass(Peptide) > ParentMass(Spectrum)
                    remove Peptide from Leaderboard
            Leaderboard ← Trim(Leaderboard, Spectrum, N)
        output LeaderPeptide
    """
    leaderboard = [""]
    leader_peptide = ""
    expander = [key for key in INTEGER_MASS.keys() if key not in ["L", "Q"]]
    while leaderboard != []:
        leaderboard = expand_peptides(leaderboard, expander=expander)
        candidates_to_remove = []
        for peptide in leaderboard:
            spec = cyclic_spectrum(peptide)
            linear_spec = linear_spectrum(peptide)
            if linear_spec[-1] == spectrum[-1]:
                if compute_score(peptide, spectrum) > compute_score(leader_peptide, spectrum):
                    leader_peptide = peptide
            elif linear_spec[-1] > spectrum[-1]:
                candidates_to_remove.append(peptide)
        leaderboard = [
            pep for pep in leaderboard if pep not in candidates_to_remove]
        leaderboard = trim(leaderboard, spectrum, N)
    return leader_peptide


def main():
    with open("/home/snopoff/Downloads/rosalind_ba4g.txt", "r") as f:
        lines = f.readlines()
    N = int(lines[0].strip())
    spectrum = list(map(int, lines[1].strip().split(" ")))
    peptide = leaderboard_cyclopeptide_sequencing(spectrum=spectrum, N=N)
    res = "-".join([str(INTEGER_MASS[acid]) for acid in peptide])
    with open("res.txt", "w") as f:
        f.write(res)


if __name__ == "__main__":
    main()
