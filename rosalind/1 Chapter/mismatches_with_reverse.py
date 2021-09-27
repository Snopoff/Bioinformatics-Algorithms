"""
#!Frequent Words with Mismatches and Reverse Complements Problem

* Extend “Find the Most Frequent Words with Mismatches in a String” to 
* find frequent words with both mismatches and reverse complements.
* @example 
    #! ACGTTGCATGTCGCATGATGCATGAGAGCT
    #! 4 1
* returns 
    #! ATGT ACAT

#? Find the most frequent k-mers (with mismatches and reverse complements) in a DNA string.

* Given: A DNA string Text as well as integers k and d.
* Return: All k-mers Pattern maximizing the sum Count_d(Text, Pattern) + Count_d(Text, Pattern_rc) over all possible k-mers.
"""


def maxMap(freqMap: dict[str, int]):
    '''
    Takes a map of strings to integers as an input and returns the maximum value of this map as output. 
    '''
    return max(freqMap.values())


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


def hamming(str1: str, str2: str):
    """
    Hamming Distance Problem: Compute the Hamming distance between two strings.

    Input: Two strings of equal length.
    Output: The Hamming distance between these strings.

    """
    assert len(str1) == len(str2)
    n = len(str1)
    res = 0
    for i in range(0, n):
        res += str1[i] != str2[i]
    return res


def neighbors(pattern: str, d: int):
    """
    Neighbors(Pattern, d)
        if d = 0
            return {Pattern}
        if |Pattern| = 1
            return {A, C, G, T}
        Neighborhood ← an empty set
        SuffixNeighbors ← Neighbors(Suffix(Pattern), d)
        for each string Text from SuffixNeighbors
            if HammingDistance(Suffix(Pattern), Text) < d
                for each nucleotide x
                    add x • Text to Neighborhood
            else
                add FirstSymbol(Pattern) • Text to Neighborhood
        return Neighborhood

    Input: A string Pattern and an integer d.
    Output: The collection of strings Neighbors(Pattern, d).
    """
    nucl = ["A", "C", "G", "T"]
    if d == 0:
        return [pattern]
    if len(pattern) == 1:
        return nucl
    neighborhood = []
    suffix = pattern[1:]
    suffixNeighborhood = neighbors(suffix, d)
    for s in suffixNeighborhood:
        if hamming(suffix, s) < d:
            for x in nucl:
                neighborhood.append(x + s)
        else:
            neighborhood.append(pattern[0] + s)

    return neighborhood


def frequent_words_with_mismatches_and_reverse_complements(text: str, k: int, d: int):
    """
    Find the most frequent k-mers (with mismatches and reverse complements) in a string.

    Input: A DNA string Text as well as integers k and d.
    Output: All k-mers Pattern maximizing the sum Count_d(Text, Pattern) + Count_d(Text, Pattern_rc) over all possible k-mers.
    """
    patterns = []
    freqMap = {}
    n = len(text)
    for i in range(0, n - k):
        pattern = text[i:i+k]
        pattern_rc = text[i:i+k]
        neighborhood = neighbors(pattern, d)
        neighborhood_rc = neighbors(pattern_rc, d)
        for j in range(len(neighborhood)):
            neighbor = neighborhood[j]
            if neighbor in freqMap:
                freqMap[neighbor] += 1
            else:
                freqMap[neighbor] = 1
        for j in range(len(neighborhood_rc)):
            neighbor_rc = neighborhood_rc[j]
            neighbor = reverse_complement(neighbor_rc)
            if neighbor in freqMap:
                freqMap[neighbor] += 1
            else:
                freqMap[neighbor] = 1
    print(freqMap)
    m = maxMap(freqMap)
    for key, value in freqMap.items():
        if value == m:
            patterns.append(key)
    return patterns


if __name__ == '__main__':
    with open("/home/snopoff/Downloads/rosalind_ba1j.txt", "r") as f:
        lines = f.readlines()
    text = lines[0].strip()
    k, d = list(map(int, lines[1].split(" ")))
    res = frequent_words_with_mismatches_and_reverse_complements(
        text, k, d)
    with open("res.txt", "w") as f:
        f.write(" ".join(res))
