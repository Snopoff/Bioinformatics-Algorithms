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
            if not last_edges:
                break
            edge = last_edges.pop(0)
            visited.update(edge)
            cycles.append([edge])
        else:
            edge = next_edge[0]
            last_edges.remove(edge)
            visited.update(edge)
            if cycles[-1][-1][-1] != edge[0]:
                cycles[-1].append(edge[::-1])
            else:
                cycles[-1].append(edge)

    for i in range(len(cycles)):
        rotation = cycles[i].index(min(cycles[i]))
        cycles[i] = cycles[i][rotation:] + cycles[i][:rotation]

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
    cycles = search_cycles(breakpoint_graph)
    return len(edges1) - len(cycles)


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
        genome[i] = np.array([-1 * element for element in genome[i]])
    return genome


def two_break_sorting(genome1: List, genome2: List):
    """
    Finds a shortest transformation of one genome into another by 2-breaks
    @param: genome1: List -- first genome
    @param: genome2: List -- second genome
        ShortestRearrangementScenario(P, Q)
            output P
            RedEdges ← ColoredEdges(P)
            BlueEdges ← ColoredEdges(Q)
            BreakpointGraph ← the graph formed by RedEdges and BlueEdges
            while BreakpointGraph has a non-trivial cycle Cycle
                (i2,i3)<-An arbitrary edge from BlueEdges in a non trivial red-blue cycle
                (i1,i2)<-An edge from RedEdges originating at node i1
                (i3,i4)<-an edge from RedEdges originating at node i3
                RedEdges ← RedEdges with edges (i1, i2) and (i3, i4) removed
                RedEdges ← RedEdges with edges (i2, i3) and (i4, i1) added
                BreakpointGraph ← the graph formed by RedEdges and BlueEdges
                P ← 2-BreakOnGenome(P, i1 , i3 , i2 , i4 )
                output P
    """
    res = [genome1]
    edges1 = colored_edges(genome1)
    edges2 = colored_edges(genome2)
    breakpoint_graph = edges1 + edges2
    cycles = search_cycles(breakpoint_graph)
    non_trivial_cycles = list(filter(lambda x: len(x) > 2, cycles))
    print(edges1)
    print(edges2)
    print()
    while non_trivial_cycles:
        cycle = non_trivial_cycles[0]
        second_edges_in_cycle = [
            edge for edge in cycle if edge in edges2 or edge[::-1] in edges2]
        i, j = second_edges_in_cycle[0]
        edge_with_i_2 = list(
            filter(lambda x: x[0] == i or x[1] == i, edges1))[0]
        edge_with_i_3 = list(
            filter(lambda x: x[0] == j or x[1] == j, edges1))[0]
        i_1, i_2 = edge_with_i_2
        i_3, i_4 = edge_with_i_3
        edges1.remove(edge_with_i_2)
        edges1.remove(edge_with_i_3)
        edges1.extend([(i_2, i_3), (i_4, i_1)])
        breakpoint_graph = edges1 + edges2
        cycles = search_cycles(breakpoint_graph)
        non_trivial_cycles = list(filter(lambda x: len(x) > 2, cycles))
        genome1 = two_break_on_genome(genome1, [i_1, i_2, i_4, i_3])
        res.append(genome1)
        print(edges1)
        print(edges2)
        print()
    return res


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


def task2():
    """
    Solve the 2-Break Sorting Problem.
    2-Break Sorting Problem: Find a shortest transformation of one genome into another by 2-breaks.

    Input: Two genomes with circular chromosomes on the same set of synteny blocks.
    Output: The sequence of genomes resulting from applying a shortest sequence of 2-breaks transforming one genome into the other.
    """
    with open("/home/snopoff/Downloads/test.txt", "r") as f:
        lines = list(map(str.strip, f.readlines()))
    genomes = [line[1:-1].split(")(") for line in lines]
    print(genomes)
    genomes = [list(map(lambda chromosome: np.array(
        list(map(int, chromosome.split(" "))), dtype=np.int32), genome)) for genome in genomes]
    sequence = two_break_sorting(*genomes)
    res = ""
    for genome in sequence:
        for chromosome in genome:
            res += "(" + print_array(chromosome) + ")"
        res += "\n"
    print(res[:-1])


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
    edges = [list(map(int, edge.split(", "))) for edge in edges]
    print(edges)
    genome = graph_to_genome(edges)
    res = ""
    for chromosome in genome:
        res += "(" + print_array(chromosome) + ")"
    with open("res.txt", "w") as f:
        f.write(res)


def task8():
    """
    Implement 2-BreakOnGenomeGraph.

    Input: The colored edges of a genome graph GenomeGraph, followed by indices i_1 , i_2 , i_3 , and i_4 .
    Output: The colored edges of the genome graph resulting from applying
            the 2-break operation 2-BreakOnGenomeGraph(GenomeGraph, i_1 , i_2 , i_3 , i_4 ).
    """
    with open("/home/snopoff/Downloads/rosalind_ba6j.txt", "r") as f:
        lines = f.readlines()

    edges = lines[0].strip()[1:-1].split("), (")
    edges = [list(map(int, edge.split(", "))) for edge in edges]
    indices = list(map(int, lines[1].strip().split(", ")))
    edges = two_break_on_genome_graph(edges, indices)
    res = ""
    for edge in edges:
        res += "(" + ', '.join(list(map(str, edge))) + "), "
    res = res[:-2]
    with open("res.txt", "w") as f:
        f.write(res)


def task9():
    """
    Implement 2-BreakOnGenome.

    Input: A genome P, followed by indices i1 , i2 , i3 , and i4 .
    Output: The genome P' resulting from applying the 2-break operation 2-BreakOnGenome(GenomeGraph i1 , i2 , i3 , i4 ).
    """
    with open("/home/snopoff/Downloads/dataset_8224_3.txt", "r") as f:
        lines = f.readlines()
    genome = lines[0].strip()[1:-1].split(")(")
    genome = list(map(lambda chromosome: np.array(
        list(map(int, chromosome.split(" "))), dtype=np.int32), genome))
    indices = list(map(int, lines[1].strip().split(", ")))
    genome_res = two_break_on_genome(genome, indices)
    res = ""
    for edge in genome_res:
        res += "(" + print_array(edge) + ")"
    with open("res.txt", "w") as f:
        f.write(res)


if __name__ == "__main__":
    task2()
