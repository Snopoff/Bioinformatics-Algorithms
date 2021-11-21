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
    cycles = []
    n = len(edges)
    active_cycle = False
    for i in range(n-1):
        edge = edges[i]
        next_edge = edges[i + 1]
        difference = abs(edge[1] - next_edge[0])
        if difference == 1:
            if active_cycle:
                cycles[-1].append(next_edge)
            else:
                cycles.append([edge, next_edge])
                active_cycle = True
        else:
            cycles.append([next_edge])

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


def search_cycles(edges: List):
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
            print("There's no next edge; forming new cycle if possible")
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
    cycles = search_cycles(breakpoint_graph)
    return len(edges1) - len(cycles)


def task1():
    """
    Solve the 2-Break Distance Problem.

    Input: Genomes P and Q.
    Output: The 2-break distance d(P, Q).
    """
    with open("/home/snopoff/Downloads/dataset_288_4.txt", "r") as f:
        lines = list(map(str.strip, f.readlines()))
    genomes = [line[1:-1].split(")(") for line in lines]
    genomes = [list(map(lambda chromosome: np.array(
        list(map(int, chromosome.split(" "))), dtype=np.int32), genome)) for genome in genomes]
    distance = two_break_distance(*genomes)
    print(distance)
    with open("res.txt", "w") as f:
        f.write(str(distance))


def task4():
    """
    Implement ChromosomeToCycle.

    Input: A chromosome Chromosome containing n synteny blocks.
    Output: The sequence Nodes of integers between 1 and 2n resulting from applying ChromosomeToCycle to Chromosome.
    """
    with open("/home/snopoff/Downloads/dataset_8222_4.txt", "r") as f:
        chromosome = np.array(
            list(map(int, f.readline().strip('()\n').split(" "))), dtype=np.int32)
    cycle = chromosome_to_cycle(chromosome)
    representation = "(" + print_array(cycle, sign=False) + ")"
    with open("res.txt", "w") as f:
        f.write(representation)


def task5():
    """
    Implement CycleToChromosome.

    Input: A sequence Nodes of integers between 1 and 2n.
    Output: The chromosome Chromosome containing n synteny blocks resulting from applying CycleToChromosome to Nodes.
    """
    with open("/home/snopoff/Downloads/dataset_8222_5.txt", "r") as f:
        cycle = np.array(
            list(map(int, f.readline().strip('()\n').split(" "))), dtype=np.int32)
    chromosome = cycle_to_chromosome(cycle)
    representation = "(" + print_array(chromosome) + ")"
    with open("res.txt", "w") as f:
        f.write(representation)


def task6():
    """
    Implement ColoredEdges.

    Input: A genome P.
    Output: The collection of colored edges in the genome graph of P in the form (x, y).
    """
    with open("/home/snopoff/Downloads/dataset_8222_7.txt", "r") as f:
        genome = f.readline().strip()[1:-1].split(")(")
        print(genome)
    genome = list(map(lambda chromosome: np.array(
        list(map(int, chromosome.split(" "))), dtype=np.int32), genome))
    edges = colored_edges(genome)
    res = ""
    for edge in edges:
        res += "(" + ', '.join(list(map(str, edge))) + "), "
    res = res[:-2]
    with open("res.txt", "w") as f:
        f.write(res)


def task7():
    """
    Implement GraphToGenome.

    Input: The colored edges ColoredEdges of a genome graph.
    Output: The genome P corresponding to this genome graph.
    """
    with open("/home/snopoff/Downloads/dataset_8222_8.txt", "r") as f:
        edges = f.readline().strip()[1:-1].split("), (")
        print(edges)
    edges = [list(map(int, edge.split(", "))) for edge in edges]
    print(edges)
    genome = graph_to_genome(edges)
    res = ""
    for chromosome in genome:
        res += "(" + print_array(chromosome) + ")"
    with open("res.txt", "w") as f:
        f.write(res)


if __name__ == "__main__":
    task1()
