"""
#! Peptide Encoding Problem
We say that a DNA string Pattern encodes an amino acid string Peptide 
if the RNA string transcribed from either Pattern or its reverse complement Pattern translates into Peptide.

Find substrings of a genome encoding a given amino acid sequence.

Given: A DNA string Text and an amino acid string Peptide.
Return: All substrings of Text encoding Peptide (if any such substrings exist).
"""

AMINO_TO_DNA_CODE = {'K': ['AAA', 'AAG'], 'N': ['AAC', 'AAT'], 'T': ['ACA', 'ACC', 'ACG', 'ACT'], 'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'], 'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'], 'I': ['ATA', 'ATC', 'ATT'], 'M': ['ATG'], 'Q': ['CAA', 'CAG'], 'H': ['CAC', 'CAT'], 'P': ['CCA', 'CCC', 'CCG', 'CCT'], 'L': [
    'CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'], 'E': ['GAA', 'GAG'], 'D': ['GAC', 'GAT'], 'A': ['GCA', 'GCC', 'GCG', 'GCT'], 'G': ['GGA', 'GGC', 'GGG', 'GGT'], 'V': ['GTA', 'GTC', 'GTG', 'GTT'], '': ['TAA', 'TAG', 'TGA'], 'Y': ['TAC', 'TAT'], 'C': ['TGC', 'TGT'], 'W': ['TGG'], 'F': ['TTC', 'TTT']}


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
    result = []
    for pos_encoding in all_possible_encodings:
        pos_encoding_rc = reverse_complement(pos_encoding)
        if pos_encoding in dna:
            result.append(pos_encoding)
        if pos_encoding_rc in dna:
            result.append(pos_encoding_rc)
    print(result)
    return result


def main():
    with open("/home/snopoff/Downloads/rosalind_ba4b.txt", "r") as f:
        lines = f.readlines()
    dna = lines[0].strip('\n')
    peptide = lines[1].strip()
    encodings = peptide_encoding(dna, peptide)
    with open("res.txt", "w") as f:
        res = "\n".join(encodings)
        f.write(res)


if __name__ == '__main__':
    main()
