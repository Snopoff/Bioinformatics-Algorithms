from typing import List
from collections import defaultdict
import numpy as np


def dp_change(money: int, coins: List[int]):
    """
    Computes minimum number of coins with denominations coins that changes money
    @param: money: int -- money needed to be changed
    @param: coints: List[int] -- denominations

      DPChange(money, Coins)
        MinNumCoins(0) ← 0
        for m ← 1 to money
            MinNumCoins(m) ← ∞
            for i ← 0 to |Coins| - 1
                if m ≥ coini
                if MinNumCoins(m - coini) + 1 < MinNumCoins(m)
                    MinNumCoins(m) ← MinNumCoins(m - coini) + 1
        output MinNumCoins(money)
    """
    min_coins = defaultdict(np.int32)
    min_coins[0] = 0
    for m in range(1, money+1):
        min_coins[m] = np.infty
        for coin in coins:
            if m >= coin:
                if min_coins[m - coin] + 1 < min_coins[m]:
                    min_coins[m] = min_coins[m - coin] + 1
    return min_coins[money]


def manhattan_tourist(n: int, m: int, down: np.array, right: np.array):
    """
    Finds the length of the longest path on given grid
    @param: n: int -- number of rows - 1
    @param: m: int -- number of columns - 1
    @param: down: np.array -- matrix of weights of vertical moves
    @param: right: np.array -- matrix of weights of horizontal moves

     ManhattanTourist(n, m, Down, Right)
        s0, 0 ← 0
        for i ← 1 to n
            si, 0 ← si-1, 0 + downi-1, 0
        for j ← 1 to m
            s0, j ← s0, j−1 + right0, j-1
        for i ← 1 to n
            for j ← 1 to m
                si, j ← max{si - 1, j + downi-1, j, si, j - 1 + righti, j-1}
        return sn, m
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


def first_task():
    """
    Code Challenge: Solve the Change Problem. The DPChange pseudocode is reproduced below for your convenience.
     Input: An integer money and an array Coins = (coin1, ..., coind).
     Output: The minimum number of coins with denominations Coins that changes money.
    """
    with open("/home/snopoff/Downloads/dataset_243_10.txt", "r") as f:
        lines = f.readlines()
    money = int(lines[0])
    coins = list(map(int, lines[1].strip().split(",")))
    min_coins_money = dp_change(money, coins)
    print(min_coins_money)
    with open("res.txt", "w") as f:
        f.write(str(min_coins_money))


def second_task():
    """
    Find the length of a longest path in the Manhattan Tourist Problem.
     Input: Integers n and m, followed by an n × (m + 1) matrix Down and an (n + 1) × m matrix Right. 
            The two matrices are separated by the "-" symbol.
     Output: The length of a longest path from source (0, 0) to sink (n, m) in the rectangular grid 
            whose edges are defined by the matrices Down and Right.
    """
    with open("/home/snopoff/Downloads/dataset_261_10.txt", "r") as f:
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


if __name__ == '__main__':
    second_task()
