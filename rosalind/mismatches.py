"""
#!Frequent Words with Mismatches Problem

* Given strings Text and Pattern as well as an integer d, we define 
*    Count_d(Text, Pattern) 
* as the total number of occurrences of Pattern in Text with at most d mismatches 
* A most frequent k-mer with up to d mismatches in Text is simply a string Pattern maximizing 
*    Count_d(Text, Pattern) 
* among all k-mers. Note that Pattern does not need to actually appear as a substring of Text
* @example 
    #! AACAAGCTGATAAACATTTAAAGAG
    #! 5 1
* returns 
    #! AAAAA

#? Find the most frequent k-mers with mismatches in a string.

* Given: A string Text as well as integers k and d.
* Return: All most frequent k-mers with up to d mismatches in Text.
"""


def maxMap(freqMap: dict[str, int]):
    '''
    Takes a map of strings to integers as an input and returns the maximum value of this map as output. 
    '''
    return max(freqMap.values())


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


def frequent_words_with_mismatches(text: str, k: int, d: int):
    """
    FrequentWordsWithMismatches(Text, k, d)
        Patterns ← an array of strings of length 0
        freqMap ← empty map
        n ← |Text|
        for i ← 0 to n - k
            Pattern ← Text(i, k)
            neighborhood ← Neighbors(Pattern, d)
            for j ← 0 to |neighborhood| - 1
                neighbor ← neighborhood[j]
                if freqMap[neighbor] doesn't exist
                    freqMap[neighbor] ← 1
                else
                    freqMap[neighbor] ← freqMap[neighbor] + 1
        m ← MaxMap(freqMap)
        for every key Pattern in freqMap
            if freqMap[Pattern] = m
                append Pattern to Patterns
        return Patterns

    Input: A string Text as well as integers k and d. (You may assume k ≤ 12 and d ≤ 3.)
    Output: All most frequent k-mers with up to d mismatches in Text.
    """
    patterns = []
    freqMap = {}
    n = len(text)
    for i in range(0, n - k):
        pattern = text[i:i+k]
        neighborhood = neighbors(pattern, d)
        for j in range(len(neighborhood)):
            neighbor = neighborhood[j]
            if neighbor in freqMap:
                freqMap[neighbor] += 1
            else:
                freqMap[neighbor] = 1
    m = maxMap(freqMap)
    for key, value in freqMap.items():
        if value == m:
            patterns.append(key)
    return patterns


if __name__ == '__main__':
    with open("/home/snopoff/Downloads/rosalind_ba1i.txt", "r") as f:
        lines = f.readlines()
    text = lines[0].strip()
    k, d = list(map(int, lines[1].split(" ")))
    res = frequent_words_with_mismatches(text, k, d)
    with open("res.txt", "w") as f:
        f.write(" ".join(res))
