"""
#! The Change Problem
Find the minimum number of coins needed to make change.

Given: An integer money and an array Coins of positive integers.
Return: The minimum number of coins with denominations Coins that changes money.
"""
from typing import List
from collections import defaultdict
import numpy as np


def dp_change(money: int, coins: List[int]):
    """
    Computes minimum number of coins with denominations coins that changes money
    @param: money: int -- money needed to be changed
    @param: coints: List[int] -- denominations
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


def main():
    with open("/home/snopoff/Downloads/rosalind_ba5a.txt", "r") as f:
        lines = f.readlines()
    money = int(lines[0])
    coins = list(map(int, lines[1].strip().split(",")))
    min_coins_money = dp_change(money, coins)
    print(min_coins_money)
    with open("res.txt", "w") as f:
        f.write(str(min_coins_money))


if __name__ == '__main__':
    main()
