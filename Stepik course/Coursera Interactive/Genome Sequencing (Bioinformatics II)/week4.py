from typing import List, Callable
from collections import defaultdict

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
        if isinstance(acid, int):
            prefix_mass[i] = prefix_mass[i-1] + acid
        else:
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
        if isinstance(acid, int):
            prefix_mass[i] = prefix_mass[i-1] + acid
        else:
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


def compute_score(peptide: str, spectrum: List[int], compute_spectrum=cyclic_spectrum):
    """
    Computes score for given peptide and spectrum
    @param: peptide: str -- given peptide
    @param: spectrum: List[str] -- given spectrum
    """
    spec = compute_spectrum(peptide)
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
    scores = sorted([compute_score(peptide, spectrum, linear_spectrum)
                    for peptide in leaderboard], reverse=False)
    if N > len(scores):
        return leaderboard
    lowest_score = scores[-N]
    result = []
    for peptide in leaderboard:
        if compute_score(peptide, spectrum, linear_spectrum) >= lowest_score:
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


def compute_convolution(spectrum: List[int]):
    """
    Computes convolution for given spectrum
    @param: spectrum: List[int] -- given spectrum 
    """
    n = len(spectrum)
    conv = []
    for i in range(0, n-1):
        for j in range(i, n):
            conv.append(spectrum[j] - spectrum[i])
    return [mass for mass in conv if mass != 0]


def choose_frequent(array: List[int], M: int, limits: List[int]):
    """
    Returns M most frequent elements of an array of integers
    @param: array: List[int] -- given array of integers
    @param: M: int -- provided number s.t. it returns M most frequent elements
    @param: limits: List[int] -- limits in which to choose elements
    """
    count = defaultdict(int)
    for element in array:
        if element >= limits[0] and element <= limits[1]:
            count[element] += 1
    multiplicities = sorted(count.values(), reverse=True)
    if M >= len(multiplicities):
        return count.keys()
    threshold = multiplicities[M]
    most_freq = []
    for key, value in count.items():
        if value >= threshold:
            most_freq.append(key)
    return most_freq


def modified_leaderboard_cyclopeptide_sequencing(spectrum: List[int], N: int, alphabet: List[int]):
    """
    Modified version of leaderboard_cyclopeptide_sequencing special for convolution_cyclopeptide_sequencing
    in order to have an opportunity to pass alphabet as parameter and work with alphabet as with list of masses
    @param: spectrum: List[int] -- given spectrum 
    @param: N: int -- number of best N scores we consider
    @param: alphabet: List[int] -- provided alphabet

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
    leaderboard = [[]]
    leader_peptide = []
    while leaderboard != []:
        leaderboard = expand_peptides(leaderboard, expander=alphabet)
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


def convolution_cyclopeptide_sequencing(spectrum: List[int], N: int, M: int):
    """
    Finds a cyclic peptide whose theoretical spectrum matches the experimental spectrum 
    with so called Leaderboard and with convolution
    @param: spectrum: List[int] -- given spectrum 
    @param: N: int -- leaderboard parameter: the number of best N scores we consider
    @param: M: int -- convolution parameter: the number of most frequent acids
    """
    conv = compute_convolution(spectrum)
    acid_masses = choose_frequent(conv, M, [57, 200])
    alphabet = [[mass] for mass in acid_masses]
    return modified_leaderboard_cyclopeptide_sequencing(spectrum, N, alphabet)


def first_task():
    """
    Cyclopeptide Scoring Problem: Compute the score of a cyclic peptide against a spectrum.

     Input: An amino acid string Peptide and a collection of integers Spectrum.
     Output: The score of Peptide against Spectrum, Score(Peptide, Spectrum).
    """
    with open("/home/snopoff/Downloads/dataset_102_3.txt", "r") as f:
        lines = f.readlines()
    peptide = lines[0].strip()
    spectrum = list(map(int, lines[1].strip().split(" ")))
    score = compute_score(peptide, spectrum)
    with open("res.txt", "w") as f:
        f.write(str(score))


def second_task():
    """
    Implement LeaderboardCyclopeptideSequencing.

     Input: An integer N and a collection of integers Spectrum.
     Output: LeaderPeptide after running LeaderboardCyclopeptideSequencing(Spectrum, N).
    """
    with open("/home/snopoff/Downloads/dataset_102_8.txt", "r") as f:
        lines = f.readlines()
    N = int(lines[0].strip())
    spectrum = list(map(int, lines[1].strip().split(" ")))
    peptide = leaderboard_cyclopeptide_sequencing(spectrum=spectrum, N=N)
    res = "-".join([str(INTEGER_MASS[acid]) for acid in peptide])
    with open("res.txt", "w") as f:
        f.write(res)


def third_task():
    """
    Spectral Convolution Problem: Compute the convolution of a spectrum.

     Input: A collection of integers Spectrum in increasing order..
     Output: The list of elements in the convolution of Spectrum. If an element has multiplicity k, 
            it should appear exactly k times; you may return the elements in any order.
    """
    with open("/home/snopoff/Downloads/dataset_104_4.txt", "r") as f:
        line = f.readline().strip()
    spectrum = list(map(int, line.split(" ")))
    conv = compute_convolution(spectrum)
    res = " ".join(list(map(str, conv)))
    with open("res.txt", "w") as f:
        f.write(res)


def fourth_task():
    """
    Implement ConvolutionCyclopeptideSequencing.

     Input: An integer M, an integer N, and a collection of (possibly repeated) integers Spectrum.
     Output: A cyclic peptide LeaderPeptide with amino acids taken only from the top M elements (and ties) of the convolution of 
            Spectrum that fall between 57 and 200, and where the size of Leaderboard is restricted to the top N (and ties).
    """
    with open("/home/snopoff/Downloads/dataset_104_7.txt", "r") as f:
        line = f.readlines()
    M = int(line[0])
    N = int(line[1])
    spectrum = sorted(list(map(int, line[2].split(" "))))
    peptide = convolution_cyclopeptide_sequencing(spectrum, N, M)
    res = "-".join([str(acid) for acid in peptide])
    with open("res.txt", "w") as f:
        f.write(res)


def fifth_task():
    """
    Compute the score of a linear peptide with respect to a spectrum.

     Input: An amino acid string Peptide and a collection of integers Spectrum.
     Output: The linear score of Peptide with respect to Spectrum, LinearScore(Peptide, Spectrum).
    """
    with open("/home/snopoff/Downloads/dataset_4913_1.txt", "r") as f:
        lines = f.readlines()
    peptide = lines[0].strip()
    spectrum = list(map(int, lines[1].split(" ")))
    lin_score = compute_score(
        peptide, spectrum, compute_spectrum=linear_spectrum)
    with open("res.txt", "w") as f:
        f.write(str(lin_score))


def sixth_task():
    """
    Implement Trim (reproduced below).

    Input: A collection of peptides Leaderboard, a collection of integers Spectrum, and an integer N.
    Output: The N highest-scoring linear peptides on Leaderboard with respect to Spectrum.
    """
    with open("/home/snopoff/Downloads/dataset_4913_3.txt", "r") as f:
        lines = f.readlines()
    leaderboard = lines[0].strip().split(" ")
    spectrum = list(map(int, lines[1].strip().split(" ")))
    N = int(lines[2])
    leaderboard_masses = [[INTEGER_MASS[acid]
                           for acid in peptide] for peptide in leaderboard]
    trimmed = trim(leaderboard_masses, spectrum, N)
    peptides = [
        leaderboard[leaderboard_masses.index(mass)] for mass in trimmed]
    with open("res.txt", "w") as f:
        f.write(" ".join(peptides))


if __name__ == "__main__":
    fourth_task()
