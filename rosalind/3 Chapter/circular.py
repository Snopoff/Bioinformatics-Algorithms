"""
#! k-Universal Circular String Problem

A k-universal circular string is a circular string that contains every possible k-mer constructed over a given alphabet.
Find a k-universal circular binary string.

Given: An integer k.
Return: A k-universal circular string. (If multiple answers exist, you may return any one.)
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


def eulerian_cycle_from_graph(graph: List[Tuple[str, List[str]]]):
    """
    Returns eulerian cycle for given graph represented as adjacency list
    @param: graph: List[Tuple[str, List[str]]] -- given adjacency list
    """
    adj_dict = from_list_to_dict(graph)
    return eulerian_cycle_from_dict(adj_dict)


def de_brujin_from_patterns(patterns: List[str]):
    """
    Construct the de Bruijn graph from a set of k-mers.
    @param: pattern: List[str] -- set of k-mers
    """
    adj_dict = defaultdict(list)
    for pattern in patterns:
        adj_dict[pattern[:-1]].append(pattern[1:])
    return list(adj_dict.items())


def create_binary_strings(k: int):
    """
    Creates binary strings
    @param: k: int -- length of binary strings
    """
    binary_strings = [None]*(2**k)
    for i in range(2**k):
        binary_strings[i] = format(i, '#0{}b'.format(k+2))[2:]
    return binary_strings


def circular_string(k: int, verbose=True):
    """
    Solves the k-Universal Circular String Problem.
    @param: k: int -- integer
    """
    binary_strings = create_binary_strings(k)
    graph = de_brujin_from_patterns(binary_strings)
    cycle = eulerian_cycle_from_graph(graph)[:-k+1]
    if verbose:
        print("Binary strings are {}".format(binary_strings))
        print("Graph is {}".format(graph))
        print("Path is {}".format(cycle))
    text = path_to_genome(cycle)
    return text


def main():
    k = int(open("/home/snopoff/Downloads/rosalind_ba3i.txt", "r").read())
    print(k)
    text = circular_string(k)
    with open("res.txt", "w") as f:
        f.write(text)


if __name__ == '__main__':
    main()
