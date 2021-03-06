from section1 import neighbors, approx_pattern_matching, hamming
import os
import numpy as np


def motif_enumeration(Dna: list, k: int, d: int):
    """
    MotifEnumeration(Dna, k, d)
        Patterns ← an empty set
        for each k-mer Pattern in Dna
            for each k-mer Pattern’ differing from Pattern by at most d mismatches
                if Pattern' appears in each string from Dna with at most d mismatches
                    add Pattern' to Patterns
        remove duplicates from Patterns
        return Patterns

    Input: Integers k and d, followed by a space-separated collection of strings Dna.
    Output: All (k, d)-motifs in Dna.
    """
    Dna = Dna.split(" ")
    print(Dna)
    patterns = []
    for dna_string in Dna:
        for i in range(0, len(dna_string) - k):
            pattern = dna_string[i:i+k]
            neighborhood = neighbors(pattern, d)
            for neighbor in neighborhood:
                is_in_all_string = True
                for dna in Dna:
                    appearances = approx_pattern_matching(neighbor, dna, d)
                    if not appearances:
                        is_in_all_string = False
                        break
                if is_in_all_string:
                    patterns.append(neighbor)

    patterns = list(set(patterns))
    return patterns


def distance_to_text(pattern: str, text: list):
    """
    Computes distance to text, where text is given as a collection of strings
    """
    res = 0
    for string in text:
        min_ham = np.infty
        for i in range(0, len(string)-len(pattern)+1):
            ham = hamming(pattern, string[i:i+len(pattern)])
            if ham < min_ham:
                min_ham = ham
        res += min_ham

    return res


def median_string(Dna: list, k: int):
    """
    MedianString(Dna, k)
        distance ← ∞
        for each k-mer Pattern from AA…AA to TT…TT
            if distance > d(Pattern, Dna)
                distance ← d(Pattern, Dna)
                Median ← Pattern
        return Median

    Input: An integer k, followed by a space-separated collection of strings Dna.
    Output: A k-mer Pattern that minimizes d(Pattern, Dna) among all possible choices of k-mers.
    """
    dist = np.infty
    median = ''
    all_possible_patterns = neighbors("A"*k, k)
    for pattern in all_possible_patterns:
        dist_to_text = distance_to_text(pattern, Dna)
        if dist_to_text < dist:
            dist = dist_to_text
            median = pattern
    return median


def compute_prob(pattern: str, profile: np.array):
    """
    Compute probability that profile generates pattern
    """
    res = 1
    letter_indices = {
        "A": 0,
        "C": 1,
        "G": 2,
        "T": 3
    }
    for index, letter in enumerate(pattern):
        res *= profile[letter_indices[letter], index]

    return res


def most_probable(text: str, k: int, profile: np.array):
    """
    Profile-most probable k-mer in Text -- a k-mer that was most likely
    to have been generated by Profile among all k-mers in Text.

    Input: A string Text, an integer k, and a 4 × k matrix Profile.
    Output: A Profile-most probable k-mer in Text.
    """
    highest_prob = 0
    res = text[0:0+k]
    for i in range(0, len(text) - k + 1):
        pattern = text[i:i+k]
        prob = compute_prob(pattern, profile)
        if prob > highest_prob:
            highest_prob = prob
            res = pattern

    return res


def generate_profile(Dna: np.array, pseudocounts=False):
    """
    Generate profile from given sequence of strings
    """
    letter_indices = {
        "A": 0,
        "C": 1,
        "G": 2,
        "T": 3
    }
    res = np.zeros((4, len(Dna[0])))
    for i in range(len(Dna)):
        for j in range(len(Dna[0])):
            res[letter_indices[Dna[i][j]], j] += 1
    if pseudocounts:
        res = res + np.ones_like(res)
        return np.divide(res, np.sum(res, axis=0))
    else:
        return res / len(Dna[0])


def compute_consensus(profile: np.array):
    """
    Compute consensus for given profile
    """
    letter_indices = {
        0: "A",
        1: "C",
        2: "G",
        3: "T"
    }
    return "".join([letter_indices[j] for j in np.argmax(profile, axis=0)])


def greedy_motif_search(Dna: list, k: int, t: int):
    """
    GreedyMotifSearch(Dna, k, t)
        BestMotifs ← motif matrix formed by first k-mers in each string from Dna
        for each k-mer Motif in the first string from Dna
            Motif1 ← Motif
            for i = 2 to t
                form Profile from motifs Motif1, …, Motifi - 1
                Motifi ← Profile-most probable k-mer in the i-th string in Dna
            Motifs ← (Motif1, …, Motift)
            if Score(Motifs) < Score(BestMotifs)
                BestMotifs ← Motifs
        return BestMotifs

    Input: Integers k and t, followed by a space-separated collection of strings Dna.
    Output: A collection of strings BestMotifs resulting from applying GreedyMotifSearch(Dna, k, t).
    """
    BestMotifs = [string[0:k]for string in Dna]
    profile_best = generate_profile(BestMotifs)
    consensus_best = compute_consensus(profile_best)
    score_best = distance_to_text(consensus_best, BestMotifs)
    for i in range(0, len(Dna[0]) - k + 1):
        Motifs = [None]*t
        Motifs[0] = Dna[0][i: i+k]
        for j in range(1, t):
            profile = generate_profile(Motifs[:j])
            print(j, Dna[j], k, profile)
            Motifs[j] = most_probable(Dna[j], k, profile)
        consensus = compute_consensus(profile)
        score = distance_to_text(consensus, Motifs)
        if score < score_best:
            BestMotifs = Motifs
            consensus_best = consensus
            score_best = score
            profile_best = profile
    return BestMotifs


def greedy_motif_search_with_pseudocounts(Dna: list, k: int, t: int):
    """
    GreedyMotifSearch(Dna, k, t)
        form a set of k-mers BestMotifs by selecting 1st k-mers in each string from Dna
        for each k-mer Motif in the first string from Dna
            Motif1 ← Motif
            for i = 2 to t
                apply Laplace's Rule of Succession to form Profile from motifs Motif1, …, Motifi-1
                Motifi ← Profile-most probable k-mer in the i-th string in Dna
            Motifs ← (Motif1, …, Motift)
            if Score(Motifs) < Score(BestMotifs)
                BestMotifs ← Motifs
        output BestMotifs

    Input: Integers k and t, followed by a space-separated collection of strings Dna.
    Output: A collection of strings BestMotifs resulting from applying GreedyMotifSearch(Dna, k, t) with pseudocounts. 
    """
    BestMotifs = [string[0:k]for string in Dna]
    profile_best = generate_profile(BestMotifs, pseudocounts=True)
    consensus_best = compute_consensus(profile_best)
    score_best = distance_to_text(consensus_best, BestMotifs)
    for i in range(0, len(Dna[0]) - k + 1):
        Motifs = [None]*t
        Motifs[0] = Dna[0][i: i+k]
        for j in range(1, t):
            profile = generate_profile(Motifs[:j], pseudocounts=True)
            #print(j, Dna[j], k, profile)
            Motifs[j] = most_probable(Dna[j], k, profile)
        consensus = compute_consensus(profile)
        score = distance_to_text(consensus, Motifs)
        if score < score_best:
            print(profile)
            BestMotifs = Motifs
            consensus_best = consensus
            score_best = score
            profile_best = profile
    return BestMotifs


def generate_motifs(profile: np.array, Dna: list):
    """
    Generate new stack of motifs using profile
    """
    Motifs = [None]*len(Dna)
    for i in range(len(Motifs)):
        Motifs[i] = most_probable(Dna[i], profile.shape[1], profile)

    return Motifs


def randomized_motif_search(Dna: list, k: int, t: int):
    """
    RandomizedMotifSearch(Dna, k, t)
        randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
        BestMotifs ← Motifs
        while forever
            Profile ← Profile(Motifs)
            Motifs ← Motifs(Profile, Dna)
            if Score(Motifs) < Score(BestMotifs)
                BestMotifs ← Motifs
            else
                return BestMotifs

    Input: Integers k and t, followed by a space-separated collection of strings Dna.
    Output: A collection BestMotifs resulting from running RandomizedMotifSearch(Dna, k, t) 1,000 times. 
            Remember to use pseudocounts!
    """
    Motifs = [None]*t
    for i in range(0, t):
        r = np.random.randint(0, len(Dna[i])-k+1)
        Motifs[i] = Dna[i][r: r+k]
    BestMotifs = Motifs
    profile_best = generate_profile(BestMotifs, pseudocounts=True)
    consensus_best = compute_consensus(profile_best)
    score_best = distance_to_text(consensus_best, BestMotifs)

    for _ in range(1000):
        profile = generate_profile(Motifs, pseudocounts=True)
        Motifs = generate_motifs(profile, Dna)
        consensus = compute_consensus(profile)
        score = distance_to_text(consensus, Motifs)
        if score < score_best:
            BestMotifs = Motifs
            consensus_best = consensus
            score_best = score
            profile_best = profile
        else:
            print(score_best)
            return BestMotifs

    print(score_best)
    return BestMotifs


def gibbs_sampler(Dna: list, k: int, t: int, N: int):
    """
    GibbsSampler(Dna, k, t, N)
        randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
        BestMotifs ← Motifs
        for j ← 1 to N
            i ← Random(t)
            Profile ← profile matrix constructed from all strings in Motifs except for Motifi
            Motifi ← Profile-randomly generated k-mer in the i-th sequence
            if Score(Motifs) < Score(BestMotifs)
                BestMotifs ← Motifs
        return BestMotifs
    """
    Motifs = [None]*t
    for i in range(0, t):
        r = np.random.randint(0, len(Dna[i])-k+1)
        Motifs[i] = Dna[i][r: r+k]
    BestMotifs = Motifs
    profile_best = generate_profile(BestMotifs, pseudocounts=True)
    consensus_best = compute_consensus(profile_best)
    score_best = distance_to_text(consensus_best, BestMotifs)

    for _ in range(N):
        i = np.random.randint(t)
        del Motifs[i]
        profile = generate_profile(Motifs, pseudocounts=True)
        motifi = most_probable(Dna[i], k, profile)
        Motifs.insert(i, motifi)
        consensus = compute_consensus(profile)
        score = distance_to_text(consensus, Motifs)
        if score < score_best:
            BestMotifs = Motifs
            consensus_best = consensus
            score_best = score
            profile_best = profile

    print(score_best)
    return BestMotifs, score_best


def first_task(curr_dir: str):
    """
    Implement MotifEnumeration (reproduced below).

    Input: Integers k and d, followed by a space-separated collection of strings Dna.
    Output: All (k, d)-motifs in Dna.
    """
    with open("/home/snopoff/Downloads/dataset_240238_8.txt", "r") as f:
        lines = f.readlines()
    str1, str2 = [line.strip() for line in lines]
    k, d = list(map(int, str1.split(" ")))

    res = motif_enumeration(str2, k, d)
    print(len(res))
    with open(curr_dir + "/res.txt", "w") as f:
        f.write(" ".join(list(map(str, res))))


def second_task(curr_dir: str):
    """
    Code Challenge: Implement MedianString.

    Input: An integer k, followed by a space-separated collection of strings Dna.
    Output: A k-mer Pattern that minimizes d(Pattern, Dna) among all possible choices of k-mers.
            (If there are multiple such strings Pattern, then you may return any one.)
    """
    with open("/home/snopoff/Downloads/dataset_240240_9 (1).txt", "r") as f:
        lines = f.readlines()
    print(lines)
    strings = [line.strip() for line in lines]
    k = int(strings[0])

    text = strings[1:]
    print("\n{}".format(text))

    res = median_string(text, k)
    print(res)
    with open(curr_dir + "/res.txt", "w") as f:
        f.write(res)


def third_task(curr_dir: str):
    """
    Profile-most Probable k-mer Problem: Find a Profile-most probable k-mer in a string.

    Input: A string Text, an integer k, and a 4 × k matrix Profile.
    Output: A Profile-most probable k-mer in Text.
    """
    with open("/home/snopoff/Downloads/dataset_240241_3.txt", "r") as f:
        lines = f.readlines()
    print(lines)
    strings = [line.strip() for line in lines]
    text = strings[0]

    k = int(strings[1])

    profile = np.array([string.split(" ")
                       for string in strings[2:]]).astype('float64')

    res = most_probable(text, k, profile)
    print(res)
    with open(curr_dir + "/res.txt", "w") as f:
        f.write(res)


#! GreedyMotifSearch
def fourth_task(curr_dir: str):
    """
    Implement GreedyMotifSearch.

    Input: Integers k and t, followed by a space-separated collection of strings Dna.
    Output: A collection of strings BestMotifs resulting from applying GreedyMotifSearch(Dna, k, t).
    If at any step you find more than one Profile-most probable k-mer in a given string, use the one occurring first.
    """
    with open("/home/snopoff/Downloads/dataset_240241_5.txt", "r") as f:
        lines = f.readlines()
    str1, str2 = [line.strip() for line in lines]
    k, t = list(map(int, str1.split(" ")))

    Dna = str2.split(" ")

    res = greedy_motif_search(Dna, k, t)
    res = " ".join(res)
    print(len(res))
    with open(curr_dir + "/res.txt", "w") as f:
        f.write("".join(list(map(str, res))))


#! GreedyMotifSearch with Pseudocounts
def fifth_task(curr_dir: str):
    """
    Implement GreedyMotifSearch with pseudocounts.

    Input: Integers k and t, followed by a space-separated collection of strings Dna.
    Output: A collection of strings BestMotifs resulting from applying GreedyMotifSearch(Dna, k, t) with pseudocounts. 
            If at any step you find more than one Profile-most probable k-mer in a given string, use the one occurring first.
    """
    with open("/home/snopoff/Downloads/dataset_240242_9.txt", "r") as f:
        lines = f.readlines()
    str1, str2 = [line.strip() for line in lines]
    k, t = list(map(int, str1.split(" ")))

    Dna = str2.split(
        " ")

    res = greedy_motif_search_with_pseudocounts(Dna, k, t)
    res = " ".join(res)
    print(res)
    with open(curr_dir + "/res.txt", "w") as f:
        f.write("".join(list(map(str, res))))


#! RandomizedMotifSearch
def sixth_task(curr_dir: str):
    """
    Implement RandomizedMotifSearch.

    Input: Integers k and t, followed by a space-separated collection of strings Dna.
    Output: A collection BestMotifs resulting from running RandomizedMotifSearch(Dna, k, t) 1,000 times. 
            Remember to use pseudocounts!
    """
    with open("/home/snopoff/Downloads/dataset_240243_5 (1).txt", "r") as f:
        lines = f.readlines()
    str1, str2 = [line.strip() for line in lines]
    k, t = list(map(int, str1.split(" ")))

    Dna = str2.split(
        " ")

    res = randomized_motif_search(Dna, k, t)
    res = " ".join(res)
    print(res)
    with open(curr_dir + "/res.txt", "w") as f:
        f.write("".join(list(map(str, res))))


#! GibbsSampler
def seventh_task(curr_dir: str):
    """
    Implement GibbsSampler.

    Input: Integers k, t, and N, followed by a space-separated collection of strings Dna.
    Output: The strings BestMotifs resulting from running GibbsSampler(Dna, k, t, N) with 20 random starts. 
            Remember to use pseudocounts!
    """
    with open("/home/snopoff/Downloads/rosalind_ba2g.txt", "r") as f:
        lines = f.readlines()
    '''
    str1, str2 = [line.strip() for line in lines]
    k, t, N = list(map(int, str1.split(" ")))

    Dna = str2.split(
        " ")
    '''

    strings = [line.strip() for line in lines]
    k, t, N = list(map(int, strings[0].split(" ")))

    Dna = strings[1:]

    best_res, best_score = "", np.infty
    for _ in range(20):
        res, score = gibbs_sampler(Dna, k, t, N)
        if score < best_score:
            best_res = res
            best_score = score
    res = " ".join(best_res)
    print(res)
    with open(curr_dir + "/res.txt", "w") as f:
        f.write("".join(list(map(str, res))))


def second_tests(*args):
    path = "/home/snopoff/Downloads/MedianString/"
    n = len([name for name in os.listdir(path + "inputs")])
    for i in range(2, n+1):
        print("-------")
        print("Stage {}".format(i))
        with open(path + "inputs/input_{}.txt".format(i), "r") as f:
            lines = f.readlines()
        str1, str2 = [line.strip() for line in lines]
        k = int(lines[0].strip())
        text = str2.split(" ")
        print("Input data is:\n k = {}\n text = {}".format(k, text))

        res = median_string(text, k)
        print(res)
        with open(path + "outputs/output_{}.txt".format(i), "r") as f:
            exp = f.read().strip()
        assert res == exp, "In {}:\n res = {}\n exp = {}".format(i, res, exp)


if __name__ == "__main__":
    curr_dir = os.getcwd()
    seventh_task(curr_dir)
