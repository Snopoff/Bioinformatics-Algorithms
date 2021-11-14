"""
#! Global Alignment Problem
Find the highest-scoring alignment between two strings using a scoring matrix.

Given: Two amino acid strings.
Return: The maximum alignment score of these strings followed by an alignment achieving this maximum score.
        Use the BLOSUM62 scoring matrix and indel penalty σ = 5.
        (If multiple alignments achieving the maximum score exist, you may return any one.)
"""
from collections import defaultdict
from typing import Dict
import sys
sys.setrecursionlimit(10000)


def global_alignment(v: str, w: str, matr: Dict[str, Dict[str, int]], sigma=5):
    """
    Solves the global alignment problem
    @param: v: str -- first string
    @param: w: str -- second string
    @param: matr: dict[str,dict[str,int]] -- scoring matrix
    @param: sigma: int -- indel penalty
    """
    n = len(v)
    m = len(w)
    s = defaultdict(lambda: (int, int))
    backtrack = defaultdict(lambda: defaultdict[int])
    s[0, 0] = 0
    backtrack[0, 0] = 0
    for i in range(1, n + 1):
        s[i, 0] = -sigma * i
    for j in range(1, m + 1):
        s[0, j] = -sigma * j
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = 0
            if v[i-1] == w[j-1]:
                match += 1
            s[i, j] = max([
                s[i-1, j] - sigma,
                s[i, j-1] - sigma,
                s[i-1, j-1] + matr[v[i-1]][w[j-1]]
            ])
            if s[i, j] == s[i-1, j] - sigma:
                backtrack[i, j] = -1
            elif s[i, j] == s[i, j-1] - sigma:
                backtrack[i, j] = 1
            elif s[i, j] == s[i-1, j-1] + matr[v[i-1]][w[j-1]]:
                backtrack[i, j] = 0
    return s[n, m], backtrack


def output_alignment(backtrack: defaultdict, v: str, i: int, j: int, first=True):
    """
    Outputs the proper alignment
    @param: backtrack: defaultdict -- given backtracks
    @param: v: str -- given first string
    @param: i: int -- start-prefix from which we output lcs
    @param: j: int -- finish-prefix to which we output lcs
    """
    if first:
        if j == 0 or i == 0:
            return ""
        if backtrack[i, j] == -1:
            return output_alignment(backtrack, v, i-1, j) + v[i-1]
        elif backtrack[i, j] == 1:
            return output_alignment(backtrack, v, i, j-1) + "-"
        else:
            return output_alignment(backtrack, v, i-1, j-1) + v[i-1]
    if not first:
        if j == 0 or i == 0:
            return ""
        if backtrack[i, j] == -1:
            return output_alignment(backtrack, v, i-1, j, False) + "-"
        elif backtrack[i, j] == 1:
            return output_alignment(backtrack, v, i, j-1, False) + v[j-1]
        else:
            return output_alignment(backtrack, v, i-1, j-1, False) + v[j-1]


def alignment(backtrack: defaultdict, str1: str, str2: str):
    """
    Returns alignment
    @param: backtrack: defaultdict -- given backtrack
    @param: str1: str -- first string
    @param: str2: str -- second string
    """
    res1 = output_alignment(backtrack, str1, len(str1), len(str2), True)
    res2 = output_alignment(backtrack, str2, len(str1), len(str2), False)
    return res1 + "\n" + res2


def main():
    with open("/home/snopoff/Downloads/rosalind_ba5e.txt", "r") as f:
        lines = f.readlines()
    str1 = lines[0].strip()
    str2 = lines[1].strip()
    blosum = {
        'A': {'A': 4, 'C': 0, 'E': -1, 'D': -2, 'G': 0, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': -2, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 0, 'W': -3, 'V': 0, 'Y': -2},
        'C': {'A': 0, 'C': 9, 'E': -4, 'D': -3, 'G': -3, 'F': -2, 'I': -1, 'H': -3, 'K': -3, 'M': -1, 'L': -1, 'N': -3, 'Q': -3, 'P': -3, 'S': -1, 'R': -3, 'T': -1, 'W': -2, 'V': -1, 'Y': -2},
        'E': {'A': -1, 'C': -4, 'E': 5, 'D': 2, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': -2, 'L': -3, 'N': 0, 'Q': 2, 'P': -1, 'S': 0, 'R': 0, 'T': -1, 'W': -3, 'V': -2, 'Y': -2},
        'D': {'A': -2, 'C': -3, 'E': 2, 'D': 6, 'G': -1, 'F': -3, 'I': -3, 'H': -1, 'K': -1, 'M': -3, 'L': -4, 'N': 1, 'Q': 0, 'P': -1, 'S': 0, 'R': -2, 'T': -1, 'W': -4, 'V': -3, 'Y': -3},
        'G': {'A': 0, 'C': -3, 'E': -2, 'D': -1, 'G': 6, 'F': -3, 'I': -4, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0, 'Q': -2, 'P': -2, 'S': 0, 'R': -2, 'T': -2, 'W': -2, 'V': -3, 'Y': -3},
        'F': {'A': -2, 'C': -2, 'E': -3, 'D': -3, 'G': -3, 'F': 6, 'I': 0, 'H': -1, 'K': -3, 'M': 0, 'L': 0, 'N': -3, 'Q': -3, 'P': -4, 'S': -2, 'R': -3, 'T': -2, 'W': 1, 'V': -1, 'Y': 3},
        'I': {'A': -1, 'C': -1, 'E': -3, 'D': -3, 'G': -4, 'F': 0, 'I': 4, 'H': -3, 'K': -3, 'M': 1, 'L': 2, 'N': -3, 'Q': -3, 'P': -3, 'S': -2, 'R': -3, 'T': -1, 'W': -3, 'V': 3, 'Y': -1},
        'H': {'A': -2, 'C': -3, 'E': 0, 'D': -1, 'G': -2, 'F': -1, 'I': -3, 'H': 8, 'K': -1, 'M': -2, 'L': -3, 'N': 1, 'Q': 0, 'P': -2, 'S': -1, 'R': 0, 'T': -2, 'W': -2, 'V': -3, 'Y': 2},
        'K': {'A': -1, 'C': -3, 'E': 1, 'D': -1, 'G': -2, 'F': -3, 'I': -3, 'H': -1, 'K': 5, 'M': -1, 'L': -2, 'N': 0, 'Q': 1, 'P': -1, 'S': 0, 'R': 2, 'T': -1, 'W': -3, 'V': -2, 'Y': -2},
        'M': {'A': -1, 'C': -1, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 1, 'H': -2, 'K': -1, 'M': 5, 'L': 2, 'N': -2, 'Q': 0, 'P': -2, 'S': -1, 'R': -1, 'T': -1, 'W': -1, 'V': 1, 'Y': -1},
        'L': {'A': -1, 'C': -1, 'E': -3, 'D': -4, 'G': -4, 'F': 0, 'I': 2, 'H': -3, 'K': -2, 'M': 2, 'L': 4, 'N': -3, 'Q': -2, 'P': -3, 'S': -2, 'R': -2, 'T': -1, 'W': -2, 'V': 1, 'Y': -1},
        'N': {'A': -2, 'C': -3, 'E': 0, 'D': 1, 'G': 0, 'F': -3, 'I': -3, 'H': 1, 'K': 0, 'M': -2, 'L': -3, 'N': 6, 'Q': 0, 'P': -2, 'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -3, 'Y': -2},
        'Q': {'A': -1, 'C': -3, 'E': 2, 'D': 0, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': 0, 'L': -2, 'N': 0, 'Q': 5, 'P': -1, 'S': 0, 'R': 1, 'T': -1, 'W': -2, 'V': -2, 'Y': -1},
        'P': {'A': -1, 'C': -3, 'E': -1, 'D': -1, 'G': -2, 'F': -4, 'I': -3, 'H': -2, 'K': -1, 'M': -2, 'L': -3, 'N': -2, 'Q': -1, 'P': 7, 'S': -1, 'R': -2, 'T': -1, 'W': -4, 'V': -2, 'Y': -3},
        'S': {'A': 1, 'C': -1, 'E': 0, 'D': 0, 'G': 0, 'F': -2, 'I': -2, 'H': -1, 'K': 0, 'M': -1, 'L': -2, 'N': 1, 'Q': 0, 'P': -1, 'S': 4, 'R': -1, 'T': 1, 'W': -3, 'V': -2, 'Y': -2},
        'R': {'A': -1, 'C': -3, 'E': 0, 'D': -2, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 2, 'M': -1, 'L': -2, 'N': 0, 'Q': 1, 'P': -2, 'S': -1, 'R': 5, 'T': -1, 'W': -3, 'V': -3, 'Y': -2},
        'T': {'A': 0, 'C': -1, 'E': -1, 'D': -1, 'G': -2, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': 0, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 5, 'W': -2, 'V': 0, 'Y': -2},
        'W': {'A': -3, 'C': -2, 'E': -3, 'D': -4, 'G': -2, 'F': 1, 'I': -3, 'H': -2, 'K': -3, 'M': -1, 'L': -2, 'N': -4, 'Q': -2, 'P': -4, 'S': -3, 'R': -3, 'T': -2, 'W': 11, 'V': -3, 'Y': 2},
        'V': {'A': 0, 'C': -1, 'E': -2, 'D': -3, 'G': -3, 'F': -1, 'I': 3, 'H': -3, 'K': -2, 'M': 1, 'L': 1, 'N': -3, 'Q': -2, 'P': -2, 'S': -2, 'R': -3, 'T': 0, 'W': -3, 'V': 4, 'Y': -1},
        'Y': {'A': -2, 'C': -2, 'E': -2, 'D': -3, 'G': -3, 'F': 3, 'I': -1, 'H': 2, 'K': -2, 'M': -1, 'L': -1, 'N': -2, 'Q': -1, 'P': -3, 'S': -2, 'R': -2, 'T': -2, 'W': 2, 'V': -1, 'Y': 7}
    }
    score, backtrack = global_alignment(str1, str2, blosum)
    align = alignment(backtrack, str1, str2)
    res = str(score) + "\n" + align
    with open("res.txt", "w") as f:
        f.write(res)


if __name__ == '__main__':
    main()
