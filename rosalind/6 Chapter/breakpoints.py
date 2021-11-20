"""
#!Number of Breakpoints Problem
Find the number of breakpoints in a permutation.

Given: A signed permutation P.
Return: The number of breakpoints in P
"""
import numpy as np


def breakpoints_number(permutation: np.array):
    """
    Counts the number of breakpoints in a permutation
    @param: permutation: np.array -- given permutation
    """
    n = permutation.shape[0]
    perm = np.zeros(n, dtype=np.int32)
    perm = np.zeros(n+2)
    perm[-1] = n + 1
    perm[1:-1] = permutation
    res = 0
    for i in range(n+1):
        if perm[i+1] - perm[i] != 1:
            res += 1
    return res


def main():
    """
    Number of Breakpoints Problem: Find the number of breakpoints in a permutation.

    Input: A permutation.
    Output: The number of breakpoints in this permutation.
    """
    with open("/home/snopoff/Downloads/rosalind_ba6b.txt", "r") as f:
        permutation = np.array(
            list(map(int, f.readline().strip()[1:-1].split(" "))), dtype=np.int32)
    num = breakpoints_number(permutation)
    with open("res.txt", "w") as f:
        f.write(str(num))


if __name__ == "__main__":
    main()
