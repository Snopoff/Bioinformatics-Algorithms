"""
#!Minimum Skew Problem

* Define the skew of a DNA string Genome, denoted Skew(Genome), as 
* the difference between the total number of occurrences of 'G' and 'C' in Genome. 
* Let Prefix_i(Genome) denote the prefix (i.e., initial substring) of Genome of length i. 
* there is an interval of Genome of length L in which Pattern appears at least t times. 
* @example 
    #! Skew(Prefixi ("CATGGGCATCGGCCATACGCC"))
* returns
    #! 0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2

#? Find a position in a genome minimizing the skew.

* Given: A DNA string Genome.
* Return: All integer(s) i minimizing Skew(Prefixi (Text)) over all values of i (from 0 to |Genome|).
"""

import numpy as np


def changeSuffix(arr: list, pos: int, value: int):
    for i in range(pos, len(arr)):
        arr[i] = value
    return arr


def calcSkew(genome: str):
    """
    If this nucleotide is G, then 
        Skewi+1(Genome) = Skewi(Genome) + 1;
    If this nucleotide is C, then 
        Skewi+1(Genome)= Skewi(Genome) – 1; 
    Otherwise, 
        Skewi+1(Genome) = Skewi(Genome).
    """
    n = len(genome)
    skew = [0]*(n+1)
    for i in range(0, n):
        if genome[i] == "G":
            changeSuffix(skew, i+1, skew[i]+1)
        elif genome[i] == "C":
            changeSuffix(skew, i+1, skew[i]-1)
    return skew


def findMinimumSkew(genome: str):
    """
    Minimum Skew Problem: Find a position in a genome where the skew diagram attains a minimum.

    Input: A DNA string Genome.
    Output: All integer(s) i minimizing Skewi (Genome) among all values of i (from 0 to |Genome|).
    """
    skew = calcSkew(genome)
    index_min = np.argmin(skew)  # the first index only
    indices = [index_min]
    for i in range(index_min+1, len(skew)):
        if skew[i] == skew[index_min]:
            indices.append(i)
    return indices


if __name__ == '__main__':
    with open("/home/snopoff/Downloads/rosalind_ba1f.txt", "r") as f:
        genome = f.read()
    res = findMinimumSkew(genome)
    with open("res.txt", "w") as f:
        f.write(" ".join(list(map(str, res))))
