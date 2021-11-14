"""
#! Local Alignment Problem
Find the highest-scoring local alignment between two strings.

Given: Two amino acid strings.
Return: The maximum score of a local alignment of the strings, followed by a local alignment of these strings 
        achieving the maximum score. Use the PAM250 scoring matrix and indel penalty σ = 5. 
        (If multiple local alignments achieving the maximum score exist, you may return any one.)
"""
from collections import defaultdict
from typing import Dict
import sys
sys.setrecursionlimit(3000)


def local_alignment(v: str, w: str, matr: Dict[str, Dict[str, int]], sigma=5):
    """
    Solves the local alignment problem
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
        s[i, 0] = 0
    for j in range(1, m + 1):
        s[0, j] = 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = 0
            if v[i-1] == w[j-1]:
                match += 1
            s[i, j] = max([
                0,
                s[i-1, j] - sigma,
                s[i, j-1] - sigma,
                s[i-1, j-1] + matr[v[i-1]][w[j-1]]
            ])
            if s[i, j] == 0:
                backtrack[i, j] = -10
            if s[i, j] == s[i-1, j] - sigma:
                backtrack[i, j] = -1
            elif s[i, j] == s[i, j-1] - sigma:
                backtrack[i, j] = 1
            elif s[i, j] == s[i-1, j-1] + matr[v[i-1]][w[j-1]]:
                backtrack[i, j] = 0
    pair, max_score = max(s.items(), key=lambda x: x[1])
    return max_score, backtrack, pair


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
        if backtrack[i, j] == -10:
            return v[j-1]
        if backtrack[i, j] == -1:
            return output_alignment(backtrack, v, i-1, j) + v[i-1]
        elif backtrack[i, j] == 1:
            return output_alignment(backtrack, v, i, j-1) + "-"
        else:
            return output_alignment(backtrack, v, i-1, j-1) + v[i-1]
    if not first:
        if j == 0 or i == 0:
            return ""
        if backtrack[i, j] == -10:
            return v[j-1]
        if backtrack[i, j] == -1:
            return output_alignment(backtrack, v, i-1, j, False) + "-"
        elif backtrack[i, j] == 1:
            return output_alignment(backtrack, v, i, j-1, False) + v[j-1]
        else:
            return output_alignment(backtrack, v, i-1, j-1, False) + v[j-1]


def alignment(backtrack: defaultdict, str1: str, str2: str, end1: int, end2: int, trim_first=False):
    """
    Returns alignment
    @param: backtrack: defaultdict -- given backtrack
    @param: str1: str -- first string
    @param: str2: str -- second string
    """
    res1 = output_alignment(backtrack, str1, end1, end2, True)
    res2 = output_alignment(backtrack, str2, end1, end2, False)
    if trim_first:
        res1 = res1[1:]
        res2 = res2[1:]
    return res1 + "\n" + res2


def main():
    with open("/home/snopoff/Downloads/rosalind_ba5f.txt", "r") as f:
        lines = f.readlines()
    str1 = lines[0].strip()
    str2 = lines[1].strip()
    pam250 = {
        'A': {'A': 2, 'C': -2, 'E': 0, 'D': 0, 'G': 1, 'F': -3, 'I': -1, 'H': -1, 'K': -1, 'M': -1, 'L': -2, 'N': 0, 'Q': 0, 'P': 1, 'S': 1, 'R': -2, 'T': 1, 'W': -6, 'V': 0, 'Y': -3},
        'C': {'A': -2, 'C': 12, 'E': -5, 'D': -5, 'G': -3, 'F': -4, 'I': -2, 'H': -3, 'K': -5, 'M': -5, 'L': -6, 'N': -4, 'Q': -5, 'P': -3, 'S': 0, 'R': -4, 'T': -2, 'W': -8, 'V': -2, 'Y': 0},
        'E': {'A': 0, 'C': -5, 'E': 4, 'D': 3, 'G': 0, 'F': -5, 'I': -2, 'H': 1, 'K': 0, 'M': -2, 'L': -3, 'N': 1, 'Q': 2, 'P': -1, 'S': 0, 'R': -1, 'T': 0, 'W': -7, 'V': -2, 'Y': -4},
        'D': {'A': 0, 'C': -5, 'E': 3, 'D': 4, 'G': 1, 'F': -6, 'I': -2, 'H': 1, 'K': 0, 'M': -3, 'L': -4, 'N': 2, 'Q': 2, 'P': -1, 'S': 0, 'R': -1, 'T': 0, 'W': -7, 'V': -2, 'Y': -4},
        'G': {'A': 1, 'C': -3, 'E': 0, 'D': 1, 'G': 5, 'F': -5, 'I': -3, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0, 'Q': -1, 'P': 0, 'S': 1, 'R': -3, 'T': 0, 'W': -7, 'V': -1, 'Y': -5},
        'F': {'A': -3, 'C': -4, 'E': -5, 'D': -6, 'G': -5, 'F': 9, 'I': 1, 'H': -2, 'K': -5, 'M': 0, 'L': 2, 'N': -3, 'Q': -5, 'P': -5, 'S': -3, 'R': -4, 'T': -3, 'W': 0, 'V': -1, 'Y': 7},
        'I': {'A': -1, 'C': -2, 'E': -2, 'D': -2, 'G': -3, 'F': 1, 'I': 5, 'H': -2, 'K': -2, 'M': 2, 'L': 2, 'N': -2, 'Q': -2, 'P': -2, 'S': -1, 'R': -2, 'T': 0, 'W': -5, 'V': 4, 'Y': -1},
        'H': {'A': -1, 'C': -3, 'E': 1, 'D': 1, 'G': -2, 'F': -2, 'I': -2, 'H': 6, 'K': 0, 'M': -2, 'L': -2, 'N': 2, 'Q': 3, 'P': 0, 'S': -1, 'R': 2, 'T': -1, 'W': -3, 'V': -2, 'Y': 0},
        'K': {'A': -1, 'C': -5, 'E': 0, 'D': 0, 'G': -2, 'F': -5, 'I': -2, 'H': 0, 'K': 5, 'M': 0, 'L': -3, 'N': 1, 'Q': 1, 'P': -1, 'S': 0, 'R': 3, 'T': 0, 'W': -3, 'V': -2, 'Y': -4},
        'M': {'A': -1, 'C': -5, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 2, 'H': -2, 'K': 0, 'M': 6, 'L': 4, 'N': -2, 'Q': -1, 'P': -2, 'S': -2, 'R': 0, 'T': -1, 'W': -4, 'V': 2, 'Y': -2},
        'L': {'A': -2, 'C': -6, 'E': -3, 'D': -4, 'G': -4, 'F': 2, 'I': 2, 'H': -2, 'K': -3, 'M': 4, 'L': 6, 'N': -3, 'Q': -2, 'P': -3, 'S': -3, 'R': -3, 'T': -2, 'W': -2, 'V': 2, 'Y': -1},
        'N': {'A': 0, 'C': -4, 'E': 1, 'D': 2, 'G': 0, 'F': -3, 'I': -2, 'H': 2, 'K': 1, 'M': -2, 'L': -3, 'N': 2, 'Q': 1, 'P': 0, 'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -2, 'Y': -2},
        'Q': {'A': 0, 'C': -5, 'E': 2, 'D': 2, 'G': -1, 'F': -5, 'I': -2, 'H': 3, 'K': 1, 'M': -1, 'L': -2, 'N': 1, 'Q': 4, 'P': 0, 'S': -1, 'R': 1, 'T': -1, 'W': -5, 'V': -2, 'Y': -4},
        'P': {'A': 1, 'C': -3, 'E': -1, 'D': -1, 'G': 0, 'F': -5, 'I': -2, 'H': 0, 'K': -1, 'M': -2, 'L': -3, 'N': 0, 'Q': 0, 'P': 6, 'S': 1, 'R': 0, 'T': 0, 'W': -6, 'V': -1, 'Y': -5},
        'S': {'A': 1, 'C': 0, 'E': 0, 'D': 0, 'G': 1, 'F': -3, 'I': -1, 'H': -1, 'K': 0, 'M': -2, 'L': -3, 'N': 1, 'Q': -1, 'P': 1, 'S': 2, 'R': 0, 'T': 1, 'W': -2, 'V': -1, 'Y': -3},
        'R': {'A': -2, 'C': -4, 'E': -1, 'D': -1, 'G': -3, 'F': -4, 'I': -2, 'H': 2, 'K': 3, 'M': 0, 'L': -3, 'N': 0, 'Q': 1, 'P': 0, 'S': 0, 'R': 6, 'T': -1, 'W': 2, 'V': -2, 'Y': -4},
        'T': {'A': 1, 'C': -2, 'E': 0, 'D': 0, 'G': 0, 'F': -3, 'I': 0, 'H': -1, 'K': 0, 'M': -1, 'L': -2, 'N': 0, 'Q': -1, 'P': 0, 'S': 1, 'R': -1, 'T': 3, 'W': -5, 'V': 0, 'Y': -3},
        'W': {'A': -6, 'C': -8, 'E': -7, 'D': -7, 'G': -7, 'F': 0, 'I': -5, 'H': -3, 'K': -3, 'M': -4, 'L': -2, 'N': -4, 'Q': -5, 'P': -6, 'S': -2, 'R': 2, 'T': -5, 'W': 17, 'V': -6, 'Y': 0},
        'V': {'A': 0, 'C': -2, 'E': -2, 'D': -2, 'G': -1, 'F': -1, 'I': 4, 'H': -2, 'K': -2, 'M': 2, 'L': 2, 'N': -2, 'Q': -2, 'P': -1, 'S': -1, 'R': -2, 'T': 0, 'W': -6, 'V': 4, 'Y': -2},
        'Y': {'A': -3, 'C': 0, 'E': -4, 'D': -4, 'G': -5, 'F': 7, 'I': -1, 'H': 0, 'K': -4, 'M': -2, 'L': -1, 'N': -2, 'Q': -4, 'P': -5, 'S': -3, 'R': -4, 'T': -3, 'W': 0, 'V': -2, 'Y': 10}
    }
    score, backtrack, pair = local_alignment(str1, str2, pam250)
    print(score)
    align = alignment(backtrack, str1, str2, pair[0], pair[1], trim_first=True)
    res = str(score) + "\n" + align
    with open("res.txt", "w") as f:
        f.write(res)


if __name__ == "__main__":
    main()
