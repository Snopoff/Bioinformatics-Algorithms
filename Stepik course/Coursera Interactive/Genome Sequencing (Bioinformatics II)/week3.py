from collections import defaultdict
from typing import DefaultDict
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
    dna_rc = reverse_complement(dna)
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
    ind = 0
    for i in range(n):
        for j in range(n - 1):
            subpeptide = cyclic_peptide[i: i+j+1]
            print(i, j, subpeptide)
            spectrum[ind] = sum([INTEGER_MASS[acid]
                                 for acid in subpeptide])
            ind += 1
    spectrum.append(sum([INTEGER_MASS[acid] for acid in peptide]))
    spectrum.append(0)
    return spectrum


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
    spectrum = sorted(calculate_spectrum(peptide))
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


if __name__ == "__main__":
    print(INTEGER_MASS)
