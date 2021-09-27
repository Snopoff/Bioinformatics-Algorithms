"""
#!Clump Finding Problem

* Given integers L and t, a string Pattern forms an (L, t)-clump inside a (larger) string Genome if 
* there is an interval of Genome of length L in which Pattern appears at least t times. 
* @example 
    #! TGCA 
* forms a (25,3)-clump in the following Genome: 
    #! gatcagcataagggtcccTGCAATGCATGACAAGCCTGCAgttgttttac.

#? Find patterns forming clumps in a string.

* Given: A string Genome, and integers k, L, and t.
* Return: All distinct k-mers forming (L, t)-clumps in Genome.
"""


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


if __name__ == "__main__":
    with open("/home/snopoff/Downloads/rosalind_ba1e.txt", "r") as f:
        lines = f.readlines()
    text = lines[0].strip()
    k, L, t = list(map(int, lines[1].strip().split(" ")))
    res = find_clumps(text, k, L, t)
    with open("res.txt", "w") as f:
        f.write(" ".join(res))
