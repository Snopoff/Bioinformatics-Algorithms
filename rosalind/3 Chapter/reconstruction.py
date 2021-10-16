"""
#! String Reconstruction Problem

Reconstruct a string from its k-mer composition

Given: An integer k followed by a list of k-mers Patterns.
Return: A string Text with k-mer composition equal to Patterns. (If multiple answers exist, you may return any one.)
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


def de_bruijn_from_patterns(patterns: List[str]):
    """
    Construct the de Bruijn graph from a set of k-mers.
    @param: pattern: List[str] -- set of k-mers
    """
    adj_dict = defaultdict(list)
    for pattern in patterns:
        adj_dict[pattern[:-1]].append(pattern[1:])
    return list(adj_dict.items())


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


def text_reconstruction(patterns: List[str], k: int, verbose=False):
    """
    Reconstructs the string by its composition
    @param: patterns: List[str] -- composition of a string
    @param: k: int -- integer s.t. patterns are k-mers
    """
    graph = de_bruijn_from_patterns(patterns)
    path = eulerian_path(graph)
    if verbose:
        print("Graph is {}".format(graph))
        print("Path is {}".format(path))
    text = path_to_genome(path)
    return text


def main():
    with open("/home/snopoff/Downloads/rosalind_ba3h.txt", "r") as f:
        lines = f.readlines()
    print(lines)
    strings = [line.strip() for line in lines]
    k = int(strings[0])
    patterns = strings[1:]
    text = text_reconstruction(patterns, k)
    with open("res.txt", "w") as f:
        f.write(text)


if __name__ == "__main__":
    main()
