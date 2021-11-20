import numpy as np


def k_sorting_reversal(k: int, permutation: np.array):
    """
    Performs reversal called k-sorting reversal, which applies to k-sorting permutation
    @param: k: int -- index k for given k-sorted permutation
    @param: permutation: np.array -- given k-sorted permutation
    """
    start = k
    finish = np.where(np.abs(permutation) == k + 1)[0][0]
    subperm = permutation[start: finish + 1]
    reverse_subperm = (-1) * subperm[::-1]
    permutation[start: finish + 1] = reverse_subperm
    return permutation


def is_element_sorted(index: int, permutation: np.array):
    """
    Returns where element is sorted or not
    @param: index: int -- index of given element of permutation
    @param: permutation: np.array -- given permutation
    """
    return index + 1 == permutation[index]


def print_array(array: np.array):
    """
    Prints array in a nice way
    @param: array: np.array -- array
    """
    res = ' '.join(('+' if i > 0 else '') + str(i) for i in array)
    return res


def greedy_sorting(permutation: np.array):
    """
    Performs greedy sorting of given permutation and returns approximated reversal distance
    @param: permutation: np.array -- given permutation

        GreedySorting(P)
            approxReversalDistance ← 0
            for k = 1 to |P|
                if element k is not sorted
                    apply the k-sorting reversal to P
                    approxReversalDistance ← approxReversalDistance + 1
                    if k-th element of P is −k
                        apply the k-sorting reversal to P
                        approxReversalDistance ← approxReversalDistance + 1
            return approxReversalDistance
    """
    dist = 0
    res = ""
    for k in range(len(permutation)):
        if not is_element_sorted(k, permutation):
            permutation = k_sorting_reversal(k, permutation)
            dist += 1
            res += print_array(permutation) + "\n"
            if k + 1 == -permutation[k]:
                permutation = k_sorting_reversal(k, permutation)
                dist += 1
                res += print_array(permutation) + "\n"
    print(permutation)
    return dist, res


def task1():
    """
    Implement GreedySorting.

    Input: A permutation P.
    Output: The sequence of permutations corresponding to applying GreedySorting to P, 
            ending with the identity permutation.
    """
    with open("/home/snopoff/Downloads/dataset_286_4.txt", "r") as f:
        permutation = np.array(
            list(map(int, f.readline().strip().split(" "))), dtype=np.int32)
    dist, res = greedy_sorting(permutation)
    with open("res.txt", "w") as f:
        f.write(res[:-1])


if __name__ == "__main__":
    task1()
