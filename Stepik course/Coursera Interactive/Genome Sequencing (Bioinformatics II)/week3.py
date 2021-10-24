from collections import defaultdict
from typing import DefaultDict, List
import tqdm


def read_codon_table(path: str):
    with open(path, "r") as f:
        lines = f.readlines()
    table = defaultdict(str)
    for line in lines:
        splitted_line = line.strip('\n').split(' ')
        table[splitted_line[0]] = splitted_line[1]
    return table


def read_integer_mass(path: str):
    with open(path, "r") as f:
        lines = f.readlines()
    table = defaultdict(str)
    for line in lines:
        splitted_line = line.strip('\n').split(' ')
        table[splitted_line[0]] = int(splitted_line[1])
    return table


GENETIC_CODE = read_codon_table(
    "/home/snopoff/Documents/MIPT/9 Semester/Bioinformatics Algorithms/data/RNA_codon_table_1.txt")

AMINO_TO_DNA_CODE = {'K': ['AAA', 'AAG'], 'N': ['AAC', 'AAT'], 'T': ['ACA', 'ACC', 'ACG', 'ACT'], 'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'], 'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'], 'I': ['ATA', 'ATC', 'ATT'], 'M': ['ATG'], 'Q': ['CAA', 'CAG'], 'H': ['CAC', 'CAT'], 'P': ['CCA', 'CCC', 'CCG', 'CCT'], 'L': [
    'CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'], 'E': ['GAA', 'GAG'], 'D': ['GAC', 'GAT'], 'A': ['GCA', 'GCC', 'GCG', 'GCT'], 'G': ['GGA', 'GGC', 'GGG', 'GGT'], 'V': ['GTA', 'GTC', 'GTG', 'GTT'], '': ['TAA', 'TAG', 'TGA'], 'Y': ['TAC', 'TAT'], 'C': ['TGC', 'TGT'], 'W': ['TGG'], 'F': ['TTC', 'TTT']}

INTEGER_MASS = read_integer_mass(
    "/home/snopoff/Documents/MIPT/9 Semester/Bioinformatics Algorithms/data/integer_mass_table.txt")


def translation(rna: str):
    """
    Performs translation
    @param: rna: str -- given rna string
    """
    peptide = ""
    for i in range(0, len(rna)-3, 3):
        codon = rna[i:i + 3]
        peptide += GENETIC_CODE[codon]
    return peptide


def transcription(dna: str):
    """
    Performs transcription
    @param: dna: str -- given dna string
    """
    return dna.replace('T', 'U')


def reverse_complement(pattern: str):
    '''
    Reverse Complement Problem: Find the reverse complement of a DNA string.
        Input: A DNA string Pattern.
        Output: The reverse complement of Pattern.
    '''
    res = ""
    for i in range(len(pattern)):
        res += (pattern[i] in 'AT')*('AT'.replace(pattern[i], "")) + \
            (pattern[i] in 'CG')*('CG'.replace(pattern[i], ""))

    return res[::-1]


def peptide_encoding(dna: str, peptide: str):
    """
    Finds substrings of a genome encoding a given amino acid sequence
    @param: dna: str -- given genome string
    @param: peptide: str -- given peptide
    """
    all_possible_encodings = AMINO_TO_DNA_CODE[peptide[0]]
    for i in range(1, len(peptide)):
        codons = AMINO_TO_DNA_CODE[peptide[i]]
        encodings = []
        for pos_encoding in all_possible_encodings:
            for codon in codons:
                encodings.append(pos_encoding + codon)
        all_possible_encodings = encodings
    print(all_possible_encodings)
    result = []
    for pos_encoding in all_possible_encodings:
        pos_encoding_rc = reverse_complement(pos_encoding)
        if pos_encoding in dna:
            result.append(pos_encoding)
        if pos_encoding_rc in dna:
            result.append(pos_encoding_rc)
    print(result)
    return result


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


def cyclic_spectrum(peptide: str, verbose=False):
    """
    Calculates cyclic spectrum for provided peptide
    This is the same as calculate_spectrum function but with different realisation, which is similar to linear_spectrum
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
    for weight in tqdm.tqdm(set(INTEGER_MASS.values())):
        amount, pos_mass_dict = bf_cyclopeptide_sequencing(
            mass - weight, pos_mass_dict)
        count += amount
    pos_mass_dict[mass] = count
    return count, pos_mass_dict


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


def mass(peptide: str):
    """
    Calculates the integer mass for given peptide
    @param: peptide: str -- given peptide
    """
    return sum([INTEGER_MASS[acid] for acid in peptide])


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
            if mass(pep) == spectrum[-1]:
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


def first_task():
    """
    Protein Translation Problem: Translate an RNA string into an amino acid string.
     Input: An RNA string Pattern and the array GeneticCode.
     Output: The translation of Pattern into an amino acid string Peptide.

    """
    with open("/home/snopoff/Downloads/dataset_96_4.txt", "r") as f:
        rna = f.readline().strip('\n')
    peptide = translation(rna)
    with open("res.txt", "w") as f:
        f.write(peptide)


def second_task():
    """
    DNA string encodes and amino acid string if 
    the RNA string transcribed from either DNA string or its rc translates to amino acid string

    Peptide Encoding Problem: Find substrings of a genome encoding a given amino acid sequence.
     Input: A DNA string Text, an amino acid string Peptide, and the array GeneticCode.
     Output: All substrings of Text encoding Peptide (if any such substrings exist).
    """
    with open("/home/snopoff/Downloads/Bacillus_brevis.txt", "r") as f:
        lines = f.readlines()
    dna = lines[0].strip('\n')
    peptide = lines[1].strip()
    encodings = peptide_encoding(dna, peptide)
    with open("res.txt", "w") as f:
        res = "\n".join(encodings)
        f.write(res)


def third_task():
    """
    Generating Theoretical Spectrum Problem: Generate the theoretical spectrum of a cyclic peptide.
     Input: An amino acid string Peptide.
     Output: Cyclospectrum(Peptide).
    """
    with open("/home/snopoff/Downloads/dataset_98_4.txt", "r") as f:
        peptide = f.readline().strip()
    spectrum = sorted(calculate_spectrum(peptide)[0])
    with open("res.txt", "w") as f:
        res = " ".join(list(map(str, spectrum)))
        f.write(res)


def fourth_task():
    """
    Counting Peptides with Given Mass Problem: Compute the number of peptides of given mass.

     Input: An integer m.
     Output: The number of linear peptides having integer mass m.
    """
    with open("/home/snopoff/Downloads/dataset_99_2.txt", "r") as f:
        mass = int(f.readline().strip())
    pos_mass_dict = defaultdict(int)
    amount, _ = bf_cyclopeptide_sequencing(
        mass, pos_mass_dict)
    with open("res.txt", "w") as f:
        f.write(str(amount))


def fifth_task():
    """
    Implement CyclopeptideSequencing (pseudocode reproduced below).

     Input: A collection of (possibly repeated) integers Spectrum corresponding to an ideal experimental spectrum.
     Output: Every amino acid string Peptide such that Cyclospectrum(Peptide) = Spectrum (if such a string exists). 
    """
    with open("/home/snopoff/Downloads/dataset_100_6 (1).txt", "r") as f:
        spectrum = [int(number) for number in f.readline().strip().split(' ')]
    peptides = cyclopeptide_sequencing(spectrum)
    masses = [[INTEGER_MASS[acid] for acid in pep] for pep in peptides]
    formatted_masses = ["-".join(list(map(str, m))) for m in masses]
    res = " ".join(formatted_masses)
    with open("res.txt", "w") as f:
        f.write(" ".join(sorted(res.split(" "))))


def sixth_task():
    """
    Implement LinearSpectrum.

     Input: An amino acid string Peptide.
     Output: The linear spectrum of Peptide.
    """
    with open("/home/snopoff/Downloads/dataset_4912_2.txt", "r") as f:
        peptide = f.readline().strip()
    linear_spec = linear_spectrum(peptide)
    with open("res.txt", "w") as f:
        res = " ".join(list(map(str, linear_spec)))
        f.write(res)


if __name__ == "__main__":
    fifth_task()
