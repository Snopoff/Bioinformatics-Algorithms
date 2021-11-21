"""
#! 2-Break Distance Problem
Find the 2-break distance between two genomes.

Given: Two genomes with circular chromosomes on the same set of synteny blocks.
Return: The 2-break distance between these two genomes.
"""
import numpy as np
from typing import List


def chromosome_to_cycle(chromosome: np.array):
    """
    Returns cycle representation of given chromosome
    @param: chromosome: np.array -- given chromosome

        ChromosomeToCycle(Chromosome)
            for j ← 0 to |Chromosome|
                i ← Chromosomej
                if i > 0
                    Nodes2j ←2i−1
                    Nodes2j+1 ← 2i
                else
                    Nodes2j ← -2i
                    Nodes2j+1 ←-2i−1
            return Nodes
    """
    n = chromosome.shape[0]
    nodes = np.zeros(2*n, dtype=np.int32)
    for j in range(n):
        i = chromosome[j]
        if i > 0:
            nodes[2*j] = 2*i - 1
            nodes[2*j + 1] = 2*i
        else:
            nodes[2*j] = -2*i
            nodes[2*j + 1] = -2*i - 1
    return nodes


def colored_edges(genome: np.array):
    """
    Returns colored edges for given genome
    @param: genome: np.array -- given genome

        ColoredEdges(P)
            Edges ← an empty set
            for each chromosome Chromosome in P
                Nodes ← ChromosomeToCycle(Chromosome)
                for j ← 1 to |Chromosome|
                    add the edge (Nodes2j-1, Nodes2j) to Edges
                add the edge (Nodes2|Chromosome|-1, Nodes0) to Edges
            return Edges
    """
    edges = []
    for chromosome in genome:
        cycle = chromosome_to_cycle(chromosome)
        n = chromosome.shape[0]
        for j in range(1, n):
            edges.append((cycle[2*j - 1], cycle[2*j]))
        edges.append((cycle[2*n - 1], cycle[0]))
    return edges


def find_cycles(edges: List):
    """
    Finds cycles in the given graph
    @param: edges: List
    """
    n = len(edges)
    cycles = [[edges[0]]]
    visited = set(edges[0])
    last_edges = edges[1:]

    while True:
        next_edge = list(
            filter(lambda edge: cycles[-1][-1][0] in edge or cycles[-1][-1][1] in edge, last_edges))
        if not next_edge:
            if not last_edges:
                break
            edge = last_edges.pop(0)
            visited.update(edge)
            cycles.append([edge])
        else:
            edge = next_edge[0]
            last_edges.remove(edge)
            visited.update(edge)
            cycles[-1].append(edge)

    return cycles


def two_break_distance(genome1: List, genome2: List):
    """
    Finds 2-break distance between genome_1 and genome_2
    @param: genome1: List -- first genome
    @param: genome2: List -- second genome
    """
    edges1 = colored_edges(genome1)
    edges2 = colored_edges(genome2)
    breakpoint_graph = edges1 + edges2
    print(breakpoint_graph)
    cycles = find_cycles(breakpoint_graph)
    return len(edges1) - len(cycles)


def main():
    with open("/home/snopoff/Downloads/rosalind_ba6c.txt", "r") as f:
        lines = list(map(str.strip, f.readlines()))
    genomes = [line[1:-1].split(")(") for line in lines]
    genomes = [list(map(lambda chromosome: np.array(
        list(map(int, chromosome.split(" "))), dtype=np.int32), genome)) for genome in genomes]
    distance = two_break_distance(*genomes)
    print(distance)
    with open("res.txt", "w") as f:
        f.write(str(distance))


if __name__ == "__main__":
    main()
