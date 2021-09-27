"""
#!Approximate Pattern Matching Problem

* We say that a k-mer Pattern appears as a substring of Text with at most d mismatches if 
* there is some k-mer substring Pattern' of Text having d or fewer mismatches with Pattern, 
*   HammingDistance(Pattern, Pattern') ≤ d. 
* @example 
    #! ATTCTGGA
    #! CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC
    #! 3
* returns 
    #! 6 7 26 27 78.

#? Find all approximate occurrences of a pattern in a string.

* Given: Strings Pattern and Text along with an integer d.
* Return: All  starting positions where Pattern appears as a substring of Text with at most d mismatches.
"""


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


def approx_pattern_matching(pattern: str, genome: str, d: int):
    """
    Approximate Pattern Matching Problem: Find all approximate occurrences of a pattern in a string.

    Input: Strings Pattern and Text along with an integer d.
    Output: All starting positions where Pattern appears as a substring of Text with at most d mismatches.
    """
    res = []
    k = len(pattern)
    for i in range(0, len(genome) - k+1):
        dist = hamming(pattern, genome[i:i+k])
        if dist <= d:
            res.append(i)

    return res


if __name__ == '__main__':
    with open("/home/snopoff/Downloads/rosalind_ba1h.txt", "r") as f:
        lines = f.readlines()
    str1, str2 = [line.strip() for line in lines[:2]]
    d = int(lines[-1].strip())

    res = approx_pattern_matching(str1, str2, d)
    with open("rosalind/res.txt", "w") as f:
        f.write(" ".join(list(map(str, res))))
