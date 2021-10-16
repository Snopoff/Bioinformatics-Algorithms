"""
#! String Reconstruction from Read-Pairs Problem

Given a string Text, a (k,d)-mer is a pair of k-mers in Text separated by distance d. 
We use the notation (Pattern_1|Pattern_2) to refer to a a (k,d)-mer whose k-mers are Pattern_1 and Pattern_2. 
The (k,d)-mer composition of Text, denoted PairedComposition_{k,d}(Text), is the collection of all (k,d)-mers in
Text (including repeated (k,d)-mers).
Reconstruct a string from its paired composition.

Given: Integers k and d followed by a collection of paired k-mers PairedReads.
Return: A string Text with (k, d)-mer composition equal to PairedReads. (If multiple answers exist, you may return any one.)
"""
from typing import List, Tuple, Dict
from collections import defaultdict


def from_list_to_dict(graph: List[Tuple[str, List[str]]]):
    """
    Returns adjacency dict from adjacency list
    @param: graph: List[Tuple[str, List[str]]] -- given adjacency list
    """
    res = defaultdict(list)
    for start, finish in graph:
        res[start] = finish
    return res


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


def eulerian_cycle_from_dict(adj_dict: Dict[str, List[str]]):
    """
    Returns eulerian cycle for given graph represented as adjacency dict
    @param: adj_dict: Dict[str, List[str]] -- given adjacency dict
    """
    res = []
    cycle = [list(adj_dict.keys())[0]]
    while cycle != []:
        next_vertices = adj_dict[cycle[-1]]
        if next_vertices == []:
            res.append(cycle.pop(-1))
        else:
            cycle.append(next_vertices.pop(-1))

    return res[::-1]


def eulerian_path(graph: List[Tuple[str, List[str]]]):
    """
    Returns eulerian path for given graph represented as adjacency list
    @param: graph: List[Tuple[str, List[str]]] -- given adjacency list
    """
    adj_dict = from_list_to_dict(graph)
    all_dict_values = [v for values in adj_dict.values() for v in values]
    start_node = graph[0][0]
    end_node = graph[0][0]
    all_vertices = list(set(all_dict_values).union(set(adj_dict.keys())))
    for key in all_vertices:
        diff = len(adj_dict[key]) - all_dict_values.count(key)
        if diff == 1:
            start_node = key
        if diff == -1:
            end_node = key
    adj_dict[end_node].append(start_node)  # balancing
    print(start_node, end_node)
    cycle = eulerian_cycle_from_dict(adj_dict)[:-1]
    if True:
        print(cycle)
        indices = [i for i, x in enumerate(cycle) if x == start_node]
        print(indices)
        indices = [i for i, x in enumerate(cycle) if x == end_node]
        print(indices)
    end = cycle.index(end_node)
    return cycle[end+1::] + cycle[:end+1:]


def paired_de_bruijn(read_pairs: List[str], verbose=True):
    """
    Constructs paired de Bruijn graph
    @param: read_pairs: List[str] -- given paired k-mers
    """
    adj_dict = defaultdict(list)
    for read in read_pairs:
        def prefix_func(x): return x[:-1]
        def suffix_func(x): return x[1:]
        prefix = "|".join(list(map(prefix_func, read.split('|'))))
        suffix = "|".join(list(map(suffix_func, read.split('|'))))
        adj_dict[prefix].append(suffix)
    if verbose:
        print(list(adj_dict.items()))
    return list(adj_dict.items())


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


def text_reconstruction_from_read_pairs(read_pairs: List[str], k: int, d: int, verbose=True):
    """
    Reconstructs text from given read_pairs
    @param: read_pairs: List[str] -- given paired k-mers
    @param: k: int -- integer s.t. we consider k-mers
    @param: d: int -- interval between paired k-mers
    """
    paired_graph = paired_de_bruijn(read_pairs)
    path = eulerian_path(paired_graph)
    genome = string_reconstruction_with_gaps(path, k, d)
    if verbose:
        # print(paired_graph)
        print("Path is {}".format(path))
        print("Genome is {}".format(genome))
    return genome


def main():
    with open("/home/snopoff/Downloads/rosalind_ba3j.txt", "r") as f:
        lines = f.readlines()
    print(lines)
    strings = [line.strip() for line in lines]
    k, d = list(map(int, strings[0].split(' ')))
    read_pairs = strings[1:]
    genome = text_reconstruction_from_read_pairs(read_pairs, k, d)
    with open("res.txt", "w") as f:
        f.write(genome)


if __name__ == '__main__':
    main()
