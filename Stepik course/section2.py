from section1 import neighbors, approx_pattern_matching
import os


def MotifEnumeration(Dna: list, k: int, d: int):
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

    res = MotifEnumeration(str2, k, d)
    print(len(res))
    with open(curr_dir + "/res.txt", "w") as f:
        f.write(" ".join(list(map(str, res))))


if __name__ == "__main__":
    curr_dir = os.getcwd()
    first_task(curr_dir)
