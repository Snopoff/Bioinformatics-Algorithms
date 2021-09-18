import numpy as np
import os


def patternCount(text: str, pattern: str):
    '''
    PatternCount(Text, Pattern)
            count ← 0
            for i ← 0 to |Text| − |Pattern|
                if Text(i, |Pattern|) = Pattern
                count ← count + 1
            return count
    '''
    count = 0

    for i in range(0, len(text) - len(pattern) + 1):
        print(text[i:i+len(pattern)], pattern)
        if text[i:i+len(pattern)] == pattern:
            count += 1

    return count


def frequentWords(text: str, k: int):
    '''
    FrequentWords(Text, k)
            FrequentPatterns ← an empty set
            for i ← 0 to |Text| − k
                Pattern ← the k-mer Text(i, k)
                Count(i) ← PatternCount(Text, Pattern)
            maxCount ← maximum value in array Count
            for i ← 0 to |Text| − k
                if Count(i) = maxCount
                    add Text(i, k) to FrequentPatterns
            remove duplicates from FrequentPatterns
            return FrequentPatterns
    '''
    frequentPatternds = []
    count = [0]*(len(text)-k)
    for i in range(0, len(text)-k):
        pattern = text[i:i+k]
        count[i] = patternCount(text, pattern)
    maxCount = max(count)
    for i in range(0, len(text)-k):
        if count(i) == maxCount:
            frequentPatternds.append(text[i:i+k])
    return list(set(frequentPatternds))


def maxMap(freqMap: dict[str, int]):
    '''
    Takes a map of strings to integers as an input and returns the maximum value of this map as output. 
    '''
    return max(freqMap.values())


def frequencyTable(text: str, k: int):
    '''
    FrequencyTable(Text, k)
        freqMap ← empty map
        n ← |Text|
        for i ← 0 to n − k
            Pattern ← Text(i, k)
            if freqMap[Pattern] doesn't exist
                freqMap[Pattern]← 1
            else
            freqMap[pattern] ←freqMap[pattern]+1 
        return freqMap
    '''
    freqMap = {}
    for i in range(0, len(text) - k):
        pattern = text[i:i+k]
        if pattern in freqMap:
            freqMap[pattern] += 1
        else:
            freqMap[pattern] = 1
    return freqMap


def betterFrequentWords(text: str, k: int):
    '''
    BetterFrequentWords(Text, k)
        FrequentPatterns ← an array of strings of length 0
        freqMap ← FrequencyTable(Text, k)
        max ← MaxMap(freqMap)
        for all strings Pattern in freqMap
            if freqMap[pattern] = max
                append Pattern to frequentPatterns
        return frequentPatterns
    '''
    frequentPatterns = []
    freqMap = frequencyTable(text, k)
    M = maxMap(freqMap)
    for key, value in freqMap.items():
        if value == M:
            frequentPatterns.append(key)

    return frequentPatterns


def reverse_complement(pattern: str):
    '''
    Reverse Complement Problem: Find the reverse complement of a DNA string.
        Input: A DNA string Pattern.
        Output: The reverse complement of Pattern.
    '''
    res = ""
    for i in range(len(pattern)):
        res += (pattern[i] in 'AT')*('AT'.replace(pattern[i], "")) + \
            (pattern[i] in 'CG')*('CG'.replace(pattern[i], ""))

    return res[::-1]


def pattern_matching(pattern: str, genome: str):
    """
    Input: Two strings, Pattern and Genome.
    Output: A collection of space-separated integers specifying all starting positions where 
            Pattern appears as a substring of Genome.
    """
    res = []
    k = len(pattern)
    for i in range(0, len(genome) - k):
        if pattern == genome[i:i+k]:
            res.append(i)

    return res


def find_clumps(text: str, k: int, L: int, t: int):
    """
    FindClumps(Text, k, L, t)
        Patterns ← an array of strings of length 0
        n ← |Text|
        for every integer i between 0 and n − L
            Window ← Text(i, L)
            freqMap ← FrequencyTable(Window, k)
            for every key s in freqMap
                if freqMap[s] ≥ t
                    append s to Patterns
        remove duplicates from Patterns
        return Patterns
    """
    patterns = []
    n = len(text)
    for i in range(0, n - L):
        window = text[i: i+L]
        freqMap = frequencyTable(window, k)
        for key, value in freqMap.items():
            if value >= t:
                patterns.append(key)
    return list(set(patterns))


def changeSuffix(arr: list, pos: int, value: int):
    for i in range(pos, len(arr)):
        arr[i] = value
    return arr


def calcSkew(genome: str):
    """
    If this nucleotide is G, then 
        Skewi+1(Genome) = Skewi(Genome) + 1;
    If this nucleotide is C, then 
        Skewi+1(Genome)= Skewi(Genome) – 1; 
    Otherwise, 
        Skewi+1(Genome) = Skewi(Genome).
    """
    n = len(genome)
    skew = [0]*(n+1)
    for i in range(0, n):
        if genome[i] == "G":
            changeSuffix(skew, i+1, skew[i]+1)
        elif genome[i] == "C":
            changeSuffix(skew, i+1, skew[i]-1)
    return skew


def findMinimumSkew(genome: str):
    """
    Minimum Skew Problem: Find a position in a genome where the skew diagram attains a minimum.

    Input: A DNA string Genome.
    Output: All integer(s) i minimizing Skewi (Genome) among all values of i (from 0 to |Genome|).
    """
    skew = calcSkew(genome)
    index_min = np.argmin(skew)  # the first index only
    indices = [index_min]
    for i in range(index_min+1, len(skew)):
        if skew[i] == skew[index_min]:
            indices.append(i)
    return indices


def hamming(str1: str, str2: str):
    """
    Hamming Distance Problem: Compute the Hamming distance between two strings.

    Input: Two strings of equal length.
    Output: The Hamming distance between these strings.

    """
    assert len(str1) == len(str2)
    n = len(str1)
    res = 0
    for i in range(0, n):
        res += str1[i] != str2[i]
    return res


def approx_pattern_matching(pattern: str, genome: str, d: int):
    """
    Approximate Pattern Matching Problem: Find all approximate occurrences of a pattern in a string.

    Input: Strings Pattern and Text along with an integer d.
    Output: All starting positions where Pattern appears as a substring of Text with at most d mismatches.
    """
    res = []
    k = len(pattern)
    for i in range(0, len(genome) - k+1):
        dist = hamming(pattern, genome[i:i+k])
        if dist <= d:
            res.append(i)

    return res


def approx_pattern_count(pattern: str, text: str, d: int):
    """
    ApproximatePatternCount(Text, Pattern, d)
        count ← 0
        for i ← 0 to |Text| − |Pattern|
            Pattern′ ← Text(i , |Pattern|)
           if HammingDistance(Pattern, Pattern′) ≤ d
                count ← count + 1
        return count

    Input: Strings Pattern and Text as well as an integer d.
    Output: Count_d(Text, Pattern).
    """
    res = 0
    k = len(pattern)
    for i in range(0, len(text) - k+1):
        dist = hamming(pattern, text[i: i+k])
        if dist <= d:
            res += 1
    return res


def neighbors(pattern: str, d: int):
    """
    Neighbors(Pattern, d)
        if d = 0
            return {Pattern}
        if |Pattern| = 1
            return {A, C, G, T}
        Neighborhood ← an empty set
        SuffixNeighbors ← Neighbors(Suffix(Pattern), d)
        for each string Text from SuffixNeighbors
            if HammingDistance(Suffix(Pattern), Text) < d
                for each nucleotide x
                    add x • Text to Neighborhood
            else
                add FirstSymbol(Pattern) • Text to Neighborhood
        return Neighborhood

    Input: A string Pattern and an integer d.
    Output: The collection of strings Neighbors(Pattern, d).
    """
    nucl = ["A", "C", "G", "T"]
    if d == 0:
        return [pattern]
    if len(pattern) == 1:
        return nucl
    neighborhood = []
    suffix = pattern[1:]
    suffixNeighborhood = neighbors(suffix, d)
    for s in suffixNeighborhood:
        if hamming(suffix, s) < d:
            for x in nucl:
                neighborhood.append(x + s)
        else:
            neighborhood.append(pattern[0] + s)

    return neighborhood


def frequent_words_with_mismatches(text: str, k: int, d: int):
    """
    FrequentWordsWithMismatches(Text, k, d)
        Patterns ← an array of strings of length 0
        freqMap ← empty map
        n ← |Text|
        for i ← 0 to n - k
            Pattern ← Text(i, k)
            neighborhood ← Neighbors(Pattern, d)
            for j ← 0 to |neighborhood| - 1
                neighbor ← neighborhood[j]
                if freqMap[neighbor] doesn't exist
                    freqMap[neighbor] ← 1
                else
                    freqMap[neighbor] ← freqMap[neighbor] + 1
        m ← MaxMap(freqMap)
        for every key Pattern in freqMap
            if freqMap[Pattern] = m
                append Pattern to Patterns
        return Patterns

    Input: A string Text as well as integers k and d. (You may assume k ≤ 12 and d ≤ 3.)
    Output: All most frequent k-mers with up to d mismatches in Text.
    """
    patterns = []
    freqMap = {}
    n = len(text)
    for i in range(0, n - k):
        pattern = text[i:i+k]
        neighborhood = neighbors(pattern, d)
        for j in range(len(neighborhood)):
            neighbor = neighborhood[j]
            if neighbor in freqMap:
                freqMap[neighbor] += 1
            else:
                freqMap[neighbor] = 1
    m = maxMap(freqMap)
    for key, value in freqMap.items():
        if value == m:
            patterns.append(key)
    return patterns


def frequent_words_with_mismatches_and_reverse_complements(text: str, k: int, d: int):
    """
    Find the most frequent k-mers (with mismatches and reverse complements) in a string.

    Input: A DNA string Text as well as integers k and d.
    Output: All k-mers Pattern maximizing the sum Count_d(Text, Pattern) + Count_d(Text, Pattern_rc) over all possible k-mers.
    """
    patterns = []
    freqMap = {}
    n = len(text)
    for i in range(0, n - k):
        pattern = text[i:i+k]
        pattern_rc = text[i:i+k]
        neighborhood = neighbors(pattern, d)
        neighborhood_rc = neighbors(pattern_rc, d)
        for j in range(len(neighborhood)):
            neighbor = neighborhood[j]
            if neighbor in freqMap:
                freqMap[neighbor] += 1
            else:
                freqMap[neighbor] = 1
        for j in range(len(neighborhood_rc)):
            neighbor_rc = neighborhood_rc[j]
            neighbor = reverse_complement(neighbor_rc)
            if neighbor in freqMap:
                freqMap[neighbor] += 1
            else:
                freqMap[neighbor] = 1
    print(freqMap)
    m = maxMap(freqMap)
    for key, value in freqMap.items():
        if value == m:
            patterns.append(key)
    return patterns


def first_task():
    '''
    Implement PatternCount (reproduced below).
        Input: Strings Text and Pattern.
        Output: Count(Text, Pattern).

            PatternCount(Text, Pattern)
            count ← 0
            for i ← 0 to |Text| − |Pattern|
                if Text(i, |Pattern|) = Pattern
                count ← count + 1
            return count
    '''

    text = input("Input text")
    pattern = input("Input pattern")
    count = patternCount(text, pattern)
    print(count)


def second_task():
    '''
    A straightforward algorithm for finding the most frequent k-mers in a string Text
        FrequentWords(Text, k)
            FrequentPatterns ← an empty set
            for i ← 0 to |Text| − k
                Pattern ← the k-mer Text(i, k)
                Count(i) ← PatternCount(Text, Pattern)
            maxCount ← maximum value in array Count
            for i ← 0 to |Text| − k
                if Count(i) = maxCount
                    add Text(i, k) to FrequentPatterns
            remove duplicates from FrequentPatterns
            return FrequentPatterns
        BetterFrequentWords(Text, k)
            FrequentPatterns ← an array of strings of length 0
            freqMap ← FrequencyTable(Text, k)
            max ← MaxMap(freqMap)
            for all strings Pattern in freqMap
                if freqMap[pattern] = max
                    append Pattern to frequentPatterns
            return frequentPatterns
    '''
    text = input("Input text:\n")
    k = int(input("Input number\n"))
    res = betterFrequentWords(text, k)
    print(" ".join(res))


def third_task(pattern: str):
    """
    Given a nucleotide p, we denote its complementary nucleotide as p*. The reverse complement of a string Pattern = p1 … pn is the string Patternrc = pn* … p1* formed by taking the complement of each nucleotide in Pattern, then reversing the resulting string. We will need the solution to the following problem throughout this chapter:

    Reverse Complement Problem: Find the reverse complement of a DNA string.

    Input: A DNA string Pattern.
    Output: Patternrc , the reverse complement of Pattern.
    """
    with open("pattern.txt", "r") as f:
        pattern = f.read().strip()
    print("------------")
    res = reverse_complement(pattern).strip(" ")
    with open("res.txt", "w") as f:
        f.write(res)
    return reverse_complement(pattern)


def fourth_task():
    """
    Code Challenge: Solve the Pattern Matching Problem.

    Input: Two strings, Pattern and Genome.
    Output: A collection of space-separated integers specifying all starting positions where 
            Pattern appears as a substring of Genome.
    """
    with open("/home/snopoff/Downloads/Vibrio_cholerae.txt", "r") as f:
        genome = f.read()
    pattern = "CTTGATCAT"
    print()
    res = pattern_matching(pattern, genome)
    with open("res.txt", "w") as f:
        f.write(" ".join(list(map(str, res))))
    print(res)


#! The Clump Finding Problem
def fifth_task():
    """
    Code Challenge: Solve the Clump Finding Problem (restated below).

    Clump Finding Problem: Find patterns forming clumps in a string.

    We defined a k-mer as a "clump" if it appears many times within a short interval of the genome. 
    More formally, given integers L and t, a k-mer Pattern forms an (L, t)-clump inside a (longer) string Genome if 
    there is an interval of Genome of length L in which this k-mer appears at least t times.

    Input: A string Genome, and integers k, L, and t.
    Output: All distinct k-mers forming (L, t)-clumps in Genome.
    """

    with open("/home/snopoff/Downloads/E_coli.txt", "r") as f:
        lines = f.readlines()
    text = lines[0].strip()
    k, L, t = (9, 500, 3)  # list(map(int, lines[1].strip().split(" ")))
    res = find_clumps(text, k, L, t)
    with open("res.txt", "w") as f:
        f.write(" ".join(res))
    '''
    text = input("Enter text:\n")
    k, L, t = list(map(int, input("Enter parameters:\n").split(" ")))
    print(find_clumps(text, k, L, t))
    '''


#! The Minimum Skew Problem
def sixth_task():
    """
    Minimum Skew Problem: Find a position in a genome where the skew diagram attains a minimum.

    Input: A DNA string Genome.
    Output: All integer(s) i minimizing Skewi (Genome) among all values of i (from 0 to |Genome|).
    """
    with open("/home/snopoff/Downloads/dataset_240220_10.txt", "r") as f:
        genome = f.read()
    #genome = "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"
    index_min = findMinimumSkew(genome)
    print(index_min)


def seventh_task():
    """
    Hamming Distance Problem: Compute the Hamming distance between two strings.

    Input: Two strings of equal length.
    Output: The Hamming distance between these strings.
    """
    with open("/home/snopoff/Downloads/dataset_240221_3.txt", "r") as f:
        lines = f.readlines()
    # [line.strip() for line in lines]
    # ["GGGCCGTTGGT", "GGACCGTTGAC"]
    str1, str2 = [line.strip() for line in lines]
    dist = hamming(str1, str2)
    print(dist)


#! Approximate Pattern Matching Problem
def eigth_task():
    """
    We say that k-mer Pattern appears as a substring of Text with at most d mismatches, i.e 
        HammingDistance(Pattern, Pattern') ≤ d

    Approximate Pattern Matching Problem: Find all approximate occurrences of a pattern in a string.

    Input: Strings Pattern and Text along with an integer d.
    Output: All starting positions where Pattern appears as a substring of Text with at most d mismatches.
    """
    with open("/home/snopoff/Downloads/dataset_240221_4.txt", "r") as f:
        lines = f.readlines()
    str1, str2 = [line.strip() for line in lines[:2]]
    d = int(lines[-1].strip())

    res = approx_pattern_matching(str1, str2, d)
    print(len(res))
    with open("Stepik course/res.txt", "w") as f:
        f.write(" ".join(list(map(str, res))))


def ninth_task():
    """
    Computing Count_d(Text, Pattern) simply requires us to compute the Hamming distance 
    between Pattern and every k-mer substring of Text, which is achieved by the following pseudocode.

       ApproximatePatternCount(Text, Pattern, d)
        count ← 0
        for i ← 0 to |Text| − |Pattern|
            Pattern′ ← Text(i , |Pattern|)
           if HammingDistance(Pattern, Pattern′) ≤ d
                count ← count + 1
        return count

    Input: Strings Pattern and Text as well as an integer d.
    Output: Count_d(Text, Pattern).
    """
    with open("/home/snopoff/Downloads/dataset_240221_6.txt", "r") as f:
        lines = f.readlines()
    str1, str2 = [line.strip() for line in lines[:2]]
    d = int(lines[-1].strip())
    res = approx_pattern_count(str1, str2, d)
    print(res)


#! Frequent Words with Mismatches Problem
def tenth_task():
    """
    Code Challenge: Solve the Frequent Words with Mismatches Problem.

    Input: A string Text as well as integers k and d. (You may assume k ≤ 12 and d ≤ 3.)
    Output: All most frequent k-mers with up to d mismatches in Text.
    """
    with open("/home/snopoff/Downloads/dataset_240221_9.txt", "r") as f:
        lines = f.readlines()
    text = lines[0].strip()
    k, d = list(map(int, lines[1].split(" ")))
    res = frequent_words_with_mismatches(text, k, d)
    with open("Stepik course/res.txt", "w") as f:
        f.write(" ".join(res))


#! Frequent Words with Mismatches and Reverse Complements Problem
def eleventh_task():
    """
    Frequent Words with Mismatches and Reverse Complements Problem: 
        Find the most frequent k-mers (with mismatches and reverse complements) in a string.

    Input: A DNA string Text as well as integers k and d.
    Output: All k-mers Pattern maximizing the sum Count_d(Text, Pattern)+ Count_d(Text, Pattern_rc) over all possible k-mers.
    """
    with open("/home/snopoff/Downloads/dataset_240221_10.txt", "r") as f:
        lines = f.readlines()
    text = lines[0].strip()
    k, d = list(map(int, lines[1].split(" ")))
    res = frequent_words_with_mismatches_and_reverse_complements(
        text, k, d)
    with open("Stepik course/res.txt", "w") as f:
        f.write(" ".join(res))


def cs_task():
    """
    Code Challenge: Implement Neighbors to find the d-neighborhood of a string.

    Input: A string Pattern and an integer d.
    Output: The collection of strings Neighbors(Pattern, d).
    """
    with open("/home/snopoff/Downloads/dataset_240229_4.txt", "r") as f:
        lines = f.readlines()
    pattern = lines[0].strip()
    d = int(lines[-1].strip())
    res = neighbors(pattern, d)
    with open("Stepik course/res.txt", "w") as f:
        f.write(" ".join(res))


def eight_tests():
    path = "/home/snopoff/Downloads/ApproximatePatternMatching/"
    n = len([name for name in os.listdir(path + "inputs")])
    for i in range(1, n+1):
        print("-------")
        print("Stage {}".format(i))
        with open(path + "inputs/input_{}.txt".format(i), "r") as f:
            lines = f.readlines()
        str1, str2 = [line.strip() for line in lines[:2]]
        d = int(lines[-1].strip())
        print("Input data is:\n str1 = {}\n str2 = {}\n d = {}".format(str1, str2, d))

        res = approx_pattern_matching(str1, str2, d)
        res = " ".join(list(map(str, res)))
        with open(path + "outputs/output_{}.txt".format(i), "r") as f:
            exp = f.read().strip()
        assert res == exp, "In {}:\n res = {}\n exp = {}".format(i, res, exp)


def tenth_tests():
    path = "/home/snopoff/Downloads/FrequentWordsMismatches/"
    n = len([name for name in os.listdir(path + "inputs")])
    for i in range(1, n+1):
        print("-------")
        print("Stage {}".format(i))
        with open(path + "inputs/input_{}.txt".format(i), "r") as f:
            lines = f.readlines()
        text = lines[0].strip()
        k, d = list(map(int, lines[-1].strip().split(" ")))
        print("Input data is:\n text = {}\n k = {}\t d = {}".format(text, k, d))

        res = frequent_words_with_mismatches(text, k, d)
        res = " ".join(list(map(str, res)))
        with open(path + "outputs/output_{}.txt".format(i), "r") as f:
            exp = f.read().strip()
        assert res == exp, "In {}:\n res = {}\n exp = {}".format(i, res, exp)


def cs_tests():
    path = "/home/snopoff/Downloads/Neighborhood/"
    n = len([name for name in os.listdir(path + "inputs")])
    for i in range(1, n+1):
        print("-------")
        print("Stage {}".format(i))
        with open(path + "inputs/input_{}.txt".format(i), "r") as f:
            lines = f.readlines()
        pattern = lines[0].strip()
        d = int(lines[-1].strip())
        print("Input data is:\n text = {}\n d = {}".format(pattern, d))

        res = neighbors(pattern, d)
        res = " ".join(list(map(str, res)))
        with open(path + "outputs/output_{}.txt".format(i), "r") as f:
            exp = f.read().strip()
        assert res == exp, "In {}:\n res = {}\n exp = {}".format(i, res, exp)


if __name__ == '__main__':
    eleventh_task()
