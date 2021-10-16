"""
#! Gapped Genome Path String Problem

Reconstruct a string from a sequence of (k,d)-mers corresponding to a path in a paired de Bruijn graph.

Given: A sequence of (k, d)-mers (a1|b1), ... , (an|bn) such that Suffix(ai|bi) = Prefix(ai+1|bi+1) for all i from 1 to n-1.
Return: A string Text where the i-th k-mer in Text is equal to Suffix(ai|bi) for all i from 1 to n, if such a string exists.
"""
from typing import List, Tuple, Dict
from collections import defaultdict


def path_to_genome(path: List[str]):
    """
    Reconstruct a string from its genome path.
    @param: path: List[str] -- A sequence path of k-mers Pattern_1, … ,Pattern_n s.t. the last k-1 symbols of Pattern_i 
                               are equal to the first k-1 symbols of Pattern_{i+1} for 1 ≤ i ≤ n-1.
    """
    res = path[0]
    for i in range(1, len(path)):
        res += path[i][-1]

    return res


def string_reconstruction_with_gaps(paired_path: List[str], k: int, d: int, verbose=True):
    """
    Reconstructs string via paired path
    @param: paired_path: List[str] -- given paired k-mers
    @param: k: int -- integer s.t. we consider k-mers
    @param: d: int -- interval between paired k-mers

    StringSpelledByGappedPatterns(GappedPatterns, k, d)
        FirstPatterns ← the sequence of initial k-mers from GappedPatterns
        SecondPatterns ← the sequence of terminal k-mers from GappedPatterns
        PrefixString ← path_to_genome(FirstPatterns, k)
        SuffixString ← path_to_genome(SecondPatterns, k)
        for i = k+d+1 to |PrefixString|
            if the i-th symbol in PrefixString does not equal the (i - k - d)-th symbol in SuffixString
                return "there is no string spelled by the gapped patterns"
        return PrefixString concatenated with the last k + d symbols of SuffixString

    """
    def first_pattern_func(x): return x.split('|')[0]
    def second_pattern_func(x): return x.split('|')[1]
    first_patterns = list(map(first_pattern_func, paired_path))
    second_patterns = list(map(second_pattern_func, paired_path))
    prefix_string = path_to_genome(first_patterns)
    suffix_string = path_to_genome(second_patterns)
    index = k + d
    overlap = prefix_string[index:] == suffix_string[:-index]
    real_size = k + d + k + len(paired_path) - 1
    if verbose:
        print(prefix_string)
        print(' '*index + suffix_string)
    if not overlap:
        return None
    return prefix_string[:index] + suffix_string


def main():
    with open("/home/snopoff/Downloads/rosalind_ba3l.txt", "r") as f:
        lines = f.readlines()
    print(lines)
    strings = [line.strip() for line in lines]
    k, d = list(map(int, strings[0].split(' ')))
    paired_path = strings[1:]
    genome = string_reconstruction_with_gaps(paired_path, k, d)
    with open("res.txt", "w") as f:
        f.write(genome)


if __name__ == '__main__':
    main()
