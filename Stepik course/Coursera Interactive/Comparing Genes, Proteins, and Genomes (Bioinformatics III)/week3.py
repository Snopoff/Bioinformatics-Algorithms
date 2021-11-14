from typing import List
import numpy as np
from operator import add


def multiple_lcs(strings: List[str]):
    """
    Finds the longest common subsequence among multiple DNA strings
    @param: strings: List[str] -- list of DNA strings
    """
    n = len(strings)
    dimensions = [len(string) + 1 for string in strings]
    score_matrix = np.zeros(dimensions, dtype=np.int32)
    backtrack = np.zeros_like(score_matrix, dtype=np.int32)

    for i in range(1, dimensions[0]):
        backtrack[i, 0, 0] = 1
    for j in range(1, dimensions[1]):
        backtrack[0, j, 0] = 2
    for k in range(1, dimensions[2]):
        backtrack[0, 0, k] = 3
    for i in range(1, dimensions[0]):
        for j in range(1, dimensions[1]):
            backtrack[i, j, 0] = 4
    for j in range(1, dimensions[1]):
        for k in range(1, dimensions[2]):
            backtrack[0, j, k] = 5
    for i in range(1, dimensions[0]):
        for k in range(1, dimensions[2]):
            backtrack[i, 0, k] = 6

    for i in range(1, dimensions[0]):
        for j in range(1, dimensions[1]):
            for k in range(1, dimensions[2]):
                match = 0
                if strings[0][i-1] == strings[1][j-1] and strings[2][k-1] == strings[1][j-1]:
                    match += 1
                score_matrix[i, j, k] = max([
                    score_matrix[i-1, j, k],
                    score_matrix[i, j-1, k],
                    score_matrix[i, j, k-1],
                    score_matrix[i-1, j-1, k],
                    score_matrix[i, j-1, k-1],
                    score_matrix[i-1, j, k-1],
                    score_matrix[i-1, j-1, k-1] + match
                ])
                if score_matrix[i, j, k] == score_matrix[i-1, j, k]:
                    backtrack[i, j, k] = 1
                if score_matrix[i, j, k] == score_matrix[i, j-1, k]:
                    backtrack[i, j, k] = 2
                if score_matrix[i, j, k] == score_matrix[i, j, k-1]:
                    backtrack[i, j, k] = 3
                if score_matrix[i, j, k] == score_matrix[i-1, j-1, k]:
                    backtrack[i, j, k] = 4
                if score_matrix[i, j, k] == score_matrix[i, j-1, k-1]:
                    backtrack[i, j, k] = 5
                if score_matrix[i, j, k] == score_matrix[i-1, j, k-1]:
                    backtrack[i, j, k] = 6
                if score_matrix[i, j, k] == score_matrix[i-1, j-1, k-1] + match:
                    backtrack[i, j, k] = 7

    return score_matrix[-1, -1, -1], backtrack


def output_lcs(strings: List[str], backtrack: np.array, endpoint: tuple):
    """
    Outputs alignment
    @param: strings: List[str] -- given list of DNA strings
    @param: backtrack: np.array -- given backtrack matrix
    @param: endpoint: tuple -- coordinate of the endpoint
    """
    if backtrack[endpoint] == 7:
        return list(map(add, output_lcs(strings, backtrack, tuple(np.subtract(endpoint, 1))), [strings[i][endpoint[i]-1] for i in range(len(strings))]))
    elif backtrack[endpoint] == 6:
        return list(map(add, output_lcs(strings, backtrack, tuple(np.subtract(endpoint, (1, 0, 1)))), [strings[0][endpoint[0]-1], '-', strings[2][endpoint[2]-1]]))
    elif backtrack[endpoint] == 5:
        return list(map(add, output_lcs(strings, backtrack, tuple(np.subtract(endpoint, (0, 1, 1)))), ['-', strings[1][endpoint[1]-1], strings[2][endpoint[2]-1]]))
    elif backtrack[endpoint] == 4:
        return list(map(add, output_lcs(strings, backtrack, tuple(np.subtract(endpoint, (1, 1, 0)))), [strings[0][endpoint[0]-1], strings[1][endpoint[1]-1], '-']))
    elif backtrack[endpoint] == 3:
        return list(map(add, output_lcs(strings, backtrack, tuple(np.subtract(endpoint, (0, 0, 1)))), ['-', '-', strings[2][endpoint[2]-1]]))
    elif backtrack[endpoint] == 2:
        return list(map(add, output_lcs(strings, backtrack, tuple(np.subtract(endpoint, (0, 1, 0)))), ['-', strings[1][endpoint[1]-1], '-']))
    elif backtrack[endpoint] == 1:
        return list(map(add, output_lcs(strings, backtrack, tuple(np.subtract(endpoint, (1, 0, 0)))), [strings[0][endpoint[0]-1], '-', '-']))
    elif backtrack[endpoint] == 0:
        return ['', '', '']


def forth_task():
    """
    Solve the Multiple Longest Common Subsequence Problem.
     Input: Three DNA strings of length at most 10.
     Output: The length of a longest common subsequence of these three strings, followed by a multiple alignment of the three strings corresponding to such an alignment.
    """
    with open("/home/snopoff/Downloads/dataset_251_5 (1).txt", "r") as f:
        lines = list(map(str.strip, f.readlines()))
    score, backtrack = multiple_lcs(lines)
    print(score)
    res = output_lcs(lines, backtrack, tuple(
        [len(string) for string in lines]))
    print(res)
    with open("res.txt", "w") as f:
        f.write("{}\n".format(int(score)) + "\n".join(res))


if __name__ == "__main__":
    forth_task()
