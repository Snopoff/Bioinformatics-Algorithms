"""
#!Length of a Longest Path in the Manhattan Tourist Problem
Find the length of a longest path in a rectangular city.

Given: Integers n and m, followed by an n × (m+1) matrix Down and an (n+1) × m matrix Right. 
       The two matrices are separated by the "-" symbol.
Return: The length of a longest path from source (0, 0) to sink (n, m) in the n × m rectangular grid 
        whose edges are defined by the matrices Down and Right.
"""
from collections import defaultdict
import numpy as np


def manhattan_tourist(n: int, m: int, down: np.array, right: np.array):
    """
    Finds the length of the longest path on given grid
    @param: n: int -- number of rows - 1
    @param: m: int -- number of columns - 1
    @param: down: np.array -- matrix of weights of vertical moves
    @param: right: np.array -- matrix of weights of horizontal moves
    """
    s = defaultdict(lambda: (int, int))
    s[0, 0] = 0
    for i in range(1, n + 1):
        s[i, 0] = s[i - 1, 0] + down[i - 1, 0]
    for j in range(1, m + 1):
        s[0, j] = s[0, j - 1] + right[0, j - 1]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s[i, j] = max(s[i - 1, j] + down[i - 1, j],
                          s[i, j - 1] + right[i, j - 1])
    return s[n, m]


def main():
    with open("/home/snopoff/Downloads/rosalind_ba5b.txt", "r") as f:
        lines = f.readlines()
    n, m = list(map(int, lines[0].strip().split(" ")))
    down = np.zeros((n, m+1), dtype=np.int32)
    right = np.zeros((n+1, m), dtype=np.int32)
    for i in range(n):
        down[i, :] = np.array(list(map(int, lines[1+i].strip().split(" "))))
    for i in range(n+1):
        right[i, :] = np.array(list(map(int, lines[n+2+i].strip().split(" "))))
    length = manhattan_tourist(n, m, down, right)
    print(length)
    with open("res.txt", "w") as f:
        f.write(str(length))


if __name__ == "__main__":
    main()
