"""
#!Compute DistanceBetweenPatternAndStrings

    * DistanceBetweenPatternAndStrings(Pattern, Dna)
    *    k ← |Pattern|
    *    distance ← 0
    *    for each string Text in Dna
    *        HammingDistance ← ∞
    *        for each k-mer Pattern’ in Text
    *            if HammingDistance > HammingDistance(Pattern, Pattern’)
    *                HammingDistance ← HammingDistance(Pattern, Pattern’)
    *        distance ← distance + HammingDistance
    *    return distance

* Given: A DNA string Pattern and a collection of DNA strings Dna.
* Return: DistanceBetweenPatternAndStrings(Pattern, Dna). 
"""
import numpy as np
import os


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


def distance_between_pattern_and_strings(pattern: str, Dna: list):
    """
    DistanceBetweenPatternAndStrings(Pattern, Dna)
        k ← |Pattern|
        distance ← 0
        for each string Text in Dna
            HammingDistance ← ∞
            for each k-mer Pattern’ in Text
                if HammingDistance > HammingDistance(Pattern, Pattern’)
                    HammingDistance ← HammingDistance(Pattern, Pattern’)
            distance ← distance + HammingDistance
        return distance
    """
    k = len(pattern)
    dist = 0
    for text in Dna:
        best_ham = np.infty
        for i in range(0, len(text) - k + 1):
            ham = hamming(pattern, text[i:i+k])
            if ham < best_ham:
                best_ham = ham
        dist += best_ham
    return dist


if __name__ == '__main__':
    with open("/home/snopoff/Downloads/rosalind_ba2h.txt", "r") as f:
        lines = f.readlines()
    strings = [line.strip() for line in lines]
    pattern = strings[0]

    Dna = strings[1].split(" ")

    res = distance_between_pattern_and_strings(pattern, Dna)
    with open(os.getcwd() + "/res.txt", "w") as f:
        f.write(str(res))
