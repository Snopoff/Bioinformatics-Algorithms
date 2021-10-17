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


def de_bruijn_from_patterns(patterns: List[str]):
    """
    Construct the de Bruijn graph from a set of k-mers.
    @param: pattern: List[str] -- set of k-mers
    """
    adj_dict = defaultdict(list)
    for pattern in patterns:
        adj_dict[pattern[:-1]].append(pattern[1:])
    return list(adj_dict.items())


def from_list_to_dict(graph: List[Tuple[str, List[str]]]):
    """
    Returns adjacency dict from adjacency list
    @param: graph: List[Tuple[str, List[str]]] -- given adjacency list
    """
    res = defaultdict(list)
    for start, finish in graph:
        res[start] = finish
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


def text_reconstruction(patterns: List[str], k: int):
    """
    Reconstructs the string by its composition
    @param: patterns: List[str] -- composition of a string
    @param: k: int -- integer s.t. patterns are k-mers
    """
    graph = de_bruijn_from_patterns(patterns)
    print("Graph is {}".format(graph))
    path = eulerian_path(graph)
    print("Path is {}".format(path))
    text = path_to_genome(path)
    return text


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
    graph = de_bruijn_from_patterns(binary_strings)
    cycle = eulerian_cycle_from_graph(graph)[:-k+1]
    if verbose:
        print("Binary strings are {}".format(binary_strings))
        print("Graph is {}".format(graph))
        print("Path is {}".format(cycle))
    text = path_to_genome(cycle)
    return text


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


def first_task():
    """
    Code Challenge: Solve the Eulerian Cycle Problem.
     Input: The adjacency list of an Eulerian directed graph.
     Output: An Eulerian cycle in this graph.
    """
    with open("/home/snopoff/Downloads/dataset_203_2 (1).txt", "r") as f:
        lines = f.readlines()
    print(lines)
    adj_strings = [line.strip().split(' -> ') for line in lines]
    adj_strings = [(string[0], string[1].split(',')) for string in adj_strings]
    cycle = eulerian_cycle_from_graph(adj_strings)
    print(cycle)
    with open("res.txt", "w") as f:
        f.write("->".join(cycle))


def second_task():
    """
    Code Challenge: Solve the Eulerian Path Problem.
     Input: The adjacency list of a directed graph that has an Eulerian path.
     Output: An Eulerian path in this graph.
    """
    with open("/home/snopoff/Downloads/dataset_203_6.txt", "r") as f:
        lines = f.readlines()
    adj_strings = [line.strip().split(' -> ') for line in lines]
    adj_strings = [(string[0], string[1].split(',')) for string in adj_strings]
    path = eulerian_path(adj_strings)
    # print(cycle)
    with open("res.txt", "w") as f:
        f.write("->".join(path))


#! String Reconstruction Problem
def third_task():
    """
    Code Challenge: Solve the String Reconstruction Problem.
     Input: An integer k followed by a list of k-mers Patterns.
     Output: A string Text with k-mer composition equal to Patterns. (If multiple answers exist, you may return any one.)
    """
    with open("/home/snopoff/Downloads/test.txt", "r") as f:
        lines = f.readlines()
    print(lines)
    strings = [line.strip() for line in lines]
    k = int(strings[0])
    patterns = strings[1:]
    text = text_reconstruction(patterns, k)
    with open("res.txt", "w") as f:
        f.write(text)


#! k-Universal Circular String Problem
def fourth_task():
    """
    Code Challenge: Solve the k-Universal Circular String Problem.
     Input: An integer k.
     Output: A k-universal circular string.
    """
    k = int(input())
    text = circular_string(k)
    with open("res.txt", "w") as f:
        f.write(text)


#! String Reconstruction from Read-Pairs Problem
def fifth_task():
    """
    Code Challenge: Solve the String Reconstruction from Read-Pairs Problem.
     Input: Integers k and d followed by a collection of paired k-mers PairedReads.
     Output: A string Text with (k, d)-mer composition equal to PairedReads.
    """
    with open("/home/snopoff/Downloads/dataset_204_16.txt", "r") as f:
        lines = f.readlines()
    print(lines)
    strings = [line.strip() for line in lines]
    k, d = list(map(int, strings[0].split(' ')))
    read_pairs = strings[1:]
    genome = text_reconstruction_from_read_pairs(read_pairs, k, d)
    with open("res.txt", "w") as f:
        f.write(genome)


#! Contig Generation Problem
def sixth_task():
    """
    Contig Generation Problem: Generate the contigs from a collection of reads (with imperfect coverage).
     Input: A collection of k-mers Patterns.
     Output: All contigs in DeBruijn(Patterns). 
    """
    with open("/home/snopoff/Downloads/dataset_205_5.txt", "r") as f:
        lines = f.readlines()
    print(lines)
    patterns = [line.strip() for line in lines]
    contigs = generate_contigs(patterns)
    with open("res.txt", "w") as f:
        f.write(" ".join(contigs))


#! Gapped Genome Path String Problem
def seventh_task():
    """
    Code Challenge: Implement StringSpelledByGappedPatterns.
     Input: Integers k and d followed by a sequence of (k,d)-mers (a_1|b_1), … , (a_n|b_n) such that Suffix(a_i|b_i) = Prefix(a_i+1|b_i+1) for 1 ≤ i ≤ n-1.
     Output: A string Text of length k+d+k+n-1 such that the i-th (k,d)-mer in Text is equal to (a_i|b_i)  for 1 ≤ i ≤ n (if such a string exists).
    """
    with open("/home/snopoff/Downloads/dataset_6206_4.txt", "r") as f:
        lines = f.readlines()
    print(lines)
    strings = [line.strip() for line in lines]
    k, d = list(map(int, strings[0].split(' ')))
    paired_path = strings[1:]
    genome = string_reconstruction_with_gaps(paired_path, k, d)
    with open("res.txt", "w") as f:
        f.write(genome)


def eigth_task():
    """
    Code Challenge: Implement MaximalNonBranchingPaths.
     Input: The adjacency list of a graph whose nodes are integers.
     Output: The collection of all maximal nonbranching paths in this graph.
    """
    with open("/home/snopoff/Downloads/dataset_6207_2.txt", "r") as f:
        lines = f.readlines()
    print(lines)
    strings = [line.strip() for line in lines]
    graph = [None]*len(strings)
    for i, string in enumerate(strings):
        splitted_string = string.split(' -> ')
        graph[i] = (splitted_string[0], splitted_string[1].split(','))
    paths = maximal_non_branching_paths(graph)
    with open('res.txt', 'w') as f:
        res = ''
        for path in paths:
            res += " -> ".join(path) + "\n"
        f.write(res[:-1])


if __name__ == "__main__":
    sixth_task()
