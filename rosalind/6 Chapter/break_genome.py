"""
#! Implement 2-BreakOnGenome
Solve the 2-Break On Genome Graph Problem.

Given: A genome P, followed by indices i, i', j, and j'.
Return: The genome P' resulting from applying the 2-break operation.
"""
import numpy as np
from typing import List


def print_array(array: np.array, sign=True):
    """
    Prints array in a nice way
    @param: array: np.array -- array
    @param: sign: bool -- if True, it prints array with sign, otherwise without
    """
    if sign:
        res = ' '.join(('+' if i > 0 else '') + str(i) for i in array)
    else:
        res = ' '.join(list(map(str, array)))
    return res


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


def cycle_to_chromosome(cycle: np.array):
    """
    Inverse to `chromosome_to_cycle` function; for given cycle, it returns corresponding chromosome
    @param: cycle: np.array -- given cycle

        CycleToChromosome(Nodes)
            for j ← 0 to |Nodes|/2
                if Nodes2j < Nodes2j+1
                    Chromosomej ← Nodes2j+1 /2
                else
                    Chromosomej ← −Nodes2j/2
            return Chromosome
    """
    n_half = cycle.shape[0] // 2
    chromosome = np.zeros(n_half, dtype=np.int32)
    for j in range(n_half):
        if cycle[2*j] < cycle[2*j + 1]:
            chromosome[j] = cycle[2*j+1] // 2
        else:
            chromosome[j] = -cycle[2*j] // 2
    return chromosome


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
    Finds cycles in given graph represented by a list of edges
    @param: edges: List -- given list of edges
    """
    cycles = [[edges[0]]]
    left_edges = edges[1:]
    while True:
        pivot = cycles[-1][-1][-1]
        desirable = pivot - 1 if pivot % 2 == 0 else pivot + 1
        if not left_edges:
            break
        next_edge = [edge for edge in left_edges if desirable in edge]
        if not next_edge:
            cycles.append([left_edges[0]])
            left_edges = left_edges[1:]
        else:
            next_edge = next_edge[0]
            if abs(cycles[-1][-1][-1] - next_edge[0]) == 1:
                cycles[-1].append(next_edge)
            else:
                cycles[-1].append(next_edge[::-1])
            left_edges.remove(next_edge)
    return cycles


def graph_to_genome(edges: List):
    """
    Inverse to `colored_edges` function, constructs genome that corresponds to given list of edges of genome graph
    @param: edges: List -- given list of edges

        GraphToGenome(GenomeGraph)
            P ← an empty set of chromosomes
            for each cycle Nodes in GenomeGraph
                Nodes ← sequence of nodes in this cycle (starting from node 1)
                Chromosome ← CycleToChromosome(Nodes)
                add Chromosome to P
            return P
    """
    genome = []
    cycles = find_cycles(edges)
    for cycle in cycles:
        nodes = np.array(
            [element for edge in cycle for element in edge], dtype=np.int32)
        nodes = np.roll(nodes, 1)
        chromosome = cycle_to_chromosome(nodes)
        genome.append(chromosome)
    return genome


def two_break_on_genome_graph(edges: List, indices: List):
    """
    Performs 2-break on GenomeGraph
    @param: edges: List -- given list of edges
    @param: indices: List -- given list of indices i_1, ..., i_4

        2-BreakOnGenomeGraph(GenomeGraph, i1 , i2 , i3 , i4)
            remove colored edges (i1, i2) and (i3, i4) from GenomeGraph
            add colored edges (i1, i3) and (i2, i4) to GenomeGraph
            return GenomeGraph
    """
    edges_to_remove = [[indices[2*i], indices[2*i+1]]
                       for i in range(len(indices) // 2)]
    for edge in edges_to_remove:
        if edge in edges:
            edges.remove(edge)
        else:
            edges.remove(edge[::-1])
    edges_to_append = list(zip(*edges_to_remove))
    edges.extend(edges_to_append)
    return edges


def two_break_on_genome(genome: np.array, indices: List):
    """
    Performs 2-break on GenomeGraph
    @param: genome: np.array -- given genome
    @param: indices: List -- given list of indices i_1, ..., i_4

        2-BreakOnGenome(P, i1 , i2 , i3 , i4 )
            GenomeGraph ← BlackEdges(P) and ColoredEdges(P)
            GenomeGraph ← 2-BreakOnGenomeGraph(GenomeGraph, i1 , i2 , i3 , i4 )
            P ← GraphToGenome(GenomeGraph)
            return P
    """
    edges = list(map(list, colored_edges(genome)))
    edges = two_break_on_genome_graph(edges, indices)
    res = graph_to_genome(edges)
    return res


def reverse_genome(genome: List):
    """
    Reverses given genome
    @param: genome: List -- given genome
    """
    n = len(genome)
    for i in range(n):
        genome[i] = genome[i][::-1]
        genome[i] = [-1 * element for element in genome[i]]
    return genome


def main():
    with open("/home/snopoff/Downloads/rosalind_ba6k (1).txt", "r") as f:
        lines = f.readlines()
    genome = lines[0].strip()[1:-1].split(")(")
    genome = list(map(lambda chromosome: np.array(
        list(map(int, chromosome.split(" "))), dtype=np.int32), genome))
    indices = list(map(int, lines[1].strip().split(", ")))
    genome_res = two_break_on_genome(genome, indices)
    genome_res = reverse_genome(genome_res)
    res = ""
    for edge in genome_res:
        res += "(" + print_array(edge) + ") "
    with open("res.txt", "w") as f:
        f.write(res)


if __name__ == "__main__":
    main()
