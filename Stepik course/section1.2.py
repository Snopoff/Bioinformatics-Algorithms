
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


if __name__ == '__main__':
    second_task()
