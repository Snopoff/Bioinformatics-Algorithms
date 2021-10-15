from typing import List


def composition(text: str, k: int):
    """
    Generates the k-mer composition of text
    @param: text: str -- given text
    @param: k: int -- interger s.t. return the function returns k-mers
    """
    amount = len(text) - k + 1
    res = [None] * amount
    for i in range(amount):
        res[i] = text[i: i + k]
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


def overlap(patterns: List[str]):
    """
    Construct the overlap graph of a collection of k-mers.
    @param: patterns: List[str] -- A collection of k-mers.
    """
    adj_list = []
    adj_dict = {}  # defaultdict(set)
    for pattern in patterns:
        if pattern[:-1] not in adj_dict:  # useless for defaultdict
            adj_dict[pattern[:-1]] = set([pattern])
        adj_dict[pattern[:-1]].add(pattern)
    for pattern in patterns:
        if pattern[1:] in adj_dict:  # useless for defaultdict
            pattern_suffixes = adj_dict[pattern[1:]]
            if pattern_suffixes:
                adj_list.append((pattern, pattern_suffixes))
    return adj_list


def de_brujin_from_patterns(patterns: List[str]):
    """
    Construct the de Bruijn graph from a set of k-mers.
    @param: pattern: List[str] -- set of k-mers
    """
    adj_dict = {}
    for pattern in patterns:
        if pattern[:-1] not in adj_dict:
            adj_dict[pattern[:-1]] = [pattern[1:]]
        else:
            adj_dict[pattern[:-1]].append(pattern[1:])
    return adj_dict.items()


def de_brujin_from_text(text: str, k: int):
    """
    Constructs the de Bruijn graph of a string.
    @param: text: str -- given text
    @param: k: int -- integer s.t. we consider k-mers
    """
    adj_dict = {}
    compositions = composition(text, k)
    # return de_brujin_from_patterns(compositions)
    for comp in compositions:
        if comp[:-1] not in adj_dict:
            adj_dict[comp[:-1]] = [comp[1:]]
        else:
            adj_dict[comp[:-1]].append(comp[1:])
    return adj_dict.items()


def first_task():
    """
    Code Challenge: Solve the String Composition Problem.
     Input: An integer k and a string Text.
     Output: Compositionk(Text) (the k-mers can be provided in any order).
    """
    with open("/home/snopoff/Downloads/dataset_197_3.txt", "r") as f:
        lines = f.readlines()
    print(lines)
    strings = [line.strip() for line in lines]
    k = int(strings[0])

    text = strings[1]
    string_comp = composition(text=text, k=k)
    with open("res.txt", "w") as f:
        f.write("\n".join(string_comp))


def second_task():
    """
    Code Challenge: Solve the String Spelled by a Genome Path Problem.
    """
    with open("/home/snopoff/Downloads/dataset_198_3.txt", "r") as f:
        lines = f.readlines()
    print(lines)
    path = [line.strip() for line in lines]

    genome = path_to_genome(path)
    with open("res.txt", "w") as f:
        f.write(genome)


def third_task():
    """
    Code Challenge: Solve the Overlap Graph Problem (restated below).
     Input: A collection Patterns of k-mers.
     Output: The overlap graph Overlap(Patterns), in the form of an adjacency list.
    """
    with open("/home/snopoff/Downloads/dataset_198_10.txt", "r") as f:
        lines = f.readlines()
    print(lines)
    patterns = [line.strip() for line in lines]

    adj_list = overlap(patterns)
    with open("res.txt", "w") as f:
        for start, finish in adj_list:
            finish_str = ",".join(finish) + "\n"
            f.write("{} -> ".format(start) + finish_str)


def fourth_task():
    """
    Code Challenge: Solve the De Bruijn Graph from a String Problem.
     Input: An integer k and a string Text.
     Output: DeBruijn_k(Text), in the form of an adjacency list.
    """
    with open("/home/snopoff/Downloads/dataset_199_6.txt", "r") as f:
        lines = f.readlines()
    print(lines)
    strings = [line.strip() for line in lines]
    k = int(strings[0])
    text = strings[1]
    adj_list = de_brujin_from_text(text=text, k=k)
    with open("res.txt", "w") as f:
        for start, finish in adj_list:
            finish_str = ",".join(finish) + "\n"
            f.write("{} -> ".format(start) + finish_str)


def fifth_task():
    """
    Code Challenge: Solve the de Bruijn Graph from k-mers Problem.
     Input: A collection of k-mers Patterns.
     Output: The adjacency list of the de Bruijn graph DeBruijn(Patterns).
    """
    with open("/home/snopoff/Downloads/dataset_200_8.txt", "r") as f:
        lines = f.readlines()
    print(lines)
    patterns = [line.strip() for line in lines]
    adj_list = de_brujin_from_patterns(patterns)
    print(adj_list)
    with open("res.txt", "w") as f:
        for start, finish in adj_list:
            finish_str = ",".join(finish) + "\n"
            f.write("{} -> ".format(start) + finish_str)


if __name__ == "__main__":
    fifth_task()
