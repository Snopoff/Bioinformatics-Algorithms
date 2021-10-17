"""
#! Contig Generation Problem

Contigs are long, contiguous segments of the genome. Biologists often settle on assembling contigs
because most assemblies still have gaps in k-mer coverage, causing the de Bruijn graph to have missing edges, 
and so the search for an Eulerian path fails.
A path in a graph is called non-branching if in(v) = out(v) = 1 for each intermediate node v of this path, i.e., 
for each node except possibly the starting and ending node of a path. 
A maximal non-branching path is a non-branching path that cannot be extended into a longer non-branching path.
Contigs correspond to strings spelled by maximal non-branching paths in the de Bruijn graph.
Generate the contigs from a collection of reads (with imperfect coverage).

Given: A collection of k-mers Patterns.
Return: All contigs in DeBruijn(Patterns). (You may return the strings in any order.)
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


def de_bruijn_from_patterns(patterns: List[str]):
    """
    Construct the de Bruijn graph from a set of k-mers.
    @param: pattern: List[str] -- set of k-mers
    """
    adj_dict = defaultdict(list)
    for pattern in patterns:
        adj_dict[pattern[:-1]].append(pattern[1:])
    return list(adj_dict.items())


def maximal_non_branching_paths(graph: List[Tuple[str, List[str]]]):
    """
    Returns maximal non branching paths for given graph
    @param: graph: List[Tuple[str, List[str]]]

    MaximalNonBranchingPaths(Graph)
        Paths ← empty list
        for each node v in Graph
            if v is not a 1-in-1-out node
                if out(v) > 0
                    for each outgoing edge (v, w) from v
                        NonBranchingPath ← the path consisting of single edge (v, w)
                        while w is a 1-in-1-out node
                            extend NonBranchingPath by the edge (w, u) 
                            w ← u
                        add NonBranchingPath to the set Paths
        for each isolated cycle Cycle in Graph
            add Cycle to Paths
        return Paths
    """
    paths = []
    adj_dict = from_list_to_dict(graph)
    print(adj_dict)
    all_dict_values = [v for values in adj_dict.values() for v in values]
    all_vertices = list(set(all_dict_values).union(set(adj_dict.keys())))
    def is_one_to_one(x): return all_dict_values.count(
        x) == 1 and len(adj_dict[x]) == 1
    all_one_to_one_nodes = list(filter(is_one_to_one, all_vertices))
    used_nodes = []
    for node in all_vertices:
        if not is_one_to_one(node):
            if len(adj_dict[node]) > 0:
                for another_node in adj_dict[node]:
                    path = [node, another_node]
                    while is_one_to_one(another_node):
                        new_node = adj_dict[another_node][0]
                        path.append(new_node)
                        another_node = new_node
                    print(path)
                    used_nodes.extend(list(set(path)))
                    paths.append(path)
    while len(set(used_nodes)) != len(set(all_vertices)):
        unused_nodes = list(set(all_vertices).difference(set(used_nodes)))
        node = unused_nodes[0]
        another_node = adj_dict[node][0]
        cycle = [node, another_node]
        while another_node != node:
            another_node = adj_dict[another_node][0]
            cycle.append(another_node)
        used_nodes.extend(list(set(cycle)))
        paths.append(cycle)
    return paths


def generate_contigs(patterns: List[str]):
    """
    Generate the contigs from a collection of reads (with imperfect coverage).
    @param: patterns: List[str] -- collection of reads
    """
    graph = de_bruijn_from_patterns(patterns)
    paths = maximal_non_branching_paths(graph)
    contigs = [None]*len(paths)
    for i, path in enumerate(paths):
        contigs[i] = path_to_genome(path)
    return contigs


def main():
    with open("/home/snopoff/Downloads/rosalind_ba3k.txt", "r") as f:
        lines = f.readlines()
    print(lines)
    patterns = [line.strip() for line in lines]
    contigs = generate_contigs(patterns)
    with open("res.txt", "w") as f:
        f.write(" ".join(contigs))


if __name__ == "__main__":
    main()
