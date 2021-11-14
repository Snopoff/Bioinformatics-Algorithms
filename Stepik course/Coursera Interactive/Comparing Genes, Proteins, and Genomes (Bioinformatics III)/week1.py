from typing import List, Dict
from collections import defaultdict
import numpy as np
import sys
sys.setrecursionlimit(1500)


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


def lcs_backtrack(v: str, w: str):
    """
    Finds the longest common subsequence of 2 strings
    @param: v: str -- first string
    @param: w: str -- second string

        LCSBackTrack(v, w)
            for i ← 0 to |v|
                si, 0 ← 0
            for j ← 0 to |w|
                s0, j ← 0
            for i ← 1 to |v|
                for j ← 1 to |w|
                    match ← 0
                    if vi-1 = wj-1
                        match ← 1
                    si, j ← max{si-1, j , si,j-1 , si-1, j-1 + match }
                    if si,j = si-1,j
                        Backtracki, j ← "↓" -1
                    else if si, j = si, j-1
                        Backtracki, j ← "→" 1
                    else if si, j = si-1, j-1 + match
                        Backtracki, j ← "↘" 0
            return Backtrack
    """
    n = len(v)
    m = len(w)
    s = defaultdict(lambda: (int, int))
    backtrack = defaultdict(lambda: defaultdict[int])
    s[0, 0] = 0
    backtrack[0, 0] = 0
    for i in range(1, n + 1):
        s[i, 0] = 0
    for j in range(1, m + 1):
        s[0, j] = 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = 0
            if v[i-1] == w[j-1]:
                match += 1
            s[i, j] = max([s[i-1, j], s[i, j-1], s[i-1, j-1] + match])
            if s[i, j] == s[i-1, j]:
                backtrack[i, j] = -1
            elif s[i, j] == s[i, j-1]:
                backtrack[i, j] = 1
            elif s[i, j] == s[i-1, j-1]:
                backtrack[i, j] = 0
    return backtrack


def output_lcs(backtrack: defaultdict, v: str, i: int, j: int):
    """
    Outputs the longest common subsequence
    @param: backtrack: defaultdict -- given backtracks
    @param: v: str -- given first string
    @param: i: int -- start-prefix from which we output lcs
    @param: j: int -- finish-prefix to which we output lcs

        OutputLCS(backtrack, v, i, j)
            if i = 0 or j = 0
                return ""
            if backtracki, j = "↓"
                return OutputLCS(backtrack, v, i - 1, j)
            else if backtracki, j = "→"
                return OutputLCS(backtrack, v, i, j - 1)
            else
                return OutputLCS(backtrack, v, i - 1, j - 1) + vi
    """
    if i == 0 or j == 0:
        return ""
    if backtrack[i, j] == -1:
        return output_lcs(backtrack, v, i-1, j)
    elif backtrack[i, j] == 1:
        return output_lcs(backtrack, v, i, j-1)
    else:
        return output_lcs(backtrack, v, i-1, j-1) + v[i-1]


def topologicalSortUtil(graph, v, visited, stack):
    visited[v] = True
    edges, _ = list(zip(*graph[v])) if graph[v] != [] else [[], []]
    for i in edges:
        if visited[i] == False:
            topologicalSortUtil(graph, i, visited, stack)
    stack.insert(0, v)


def topological_sort(graph, length):
    visited = [False]*length
    stack = []
    for i in range(length):
        if visited[i] == False:
            topologicalSortUtil(graph, i, visited, stack)
    return stack


def longest_path_in_dag(graph: defaultdict[list], source: int, sink: int):
    """
    Finds the longest path in DAG
    @param: graph: defaultdict[list] -- given DAG
    @param: source: int -- start node
    @param: sink: int -- finish node

        LongestPath(Graph, source, sink)
            for each node b in Graph
                s_b ← −∞
            s_source ← 0
            topologically order Graph
            for each node b in Graph (following the topological order)
                s_b ← max_{all predecessors a of node b} {s_a + weight of edge from a to b}
            return s_sink
    """
    s = defaultdict(np.int32)
    backtrack = defaultdict(int)
    for key in graph.keys():
        s[key] = -np.infty
    s[source] = 0
    order = topological_sort(graph, sink+1)
    for vertex in order:
        edges = graph[vertex]
        for edge, weight in edges:
            value = s[vertex] + weight
            if s[edge] < value:
                s[edge] = value
                backtrack[edge] = vertex
    return s[sink], backtrack


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


def third_task():
    """
    Use OutputLCS (reproduced below) to solve the Longest Common Subsequence Problem.
     Input: Two strings s and t.
     Output: A longest common subsequence of s and t. (Note: more than one solution may exist, in which case you may output any one.)
    """
    with open("/home/snopoff/Downloads/dataset_245_5.txt", "r") as f:
        lines = f.readlines()
    str1 = lines[0].strip()
    str2 = lines[1].strip()
    backtrack = lcs_backtrack(str1, str2)
    lcs = output_lcs(backtrack, str1, len(str1), len(str2))
    with open("res.txt", "w") as f:
        f.write(lcs)


def forth_task():
    """
    Solve the Longest Path in a DAG Problem.
     Input: An integer representing the starting node to consider in a graph, followed by an integer representing 
            the ending node to consider, followed by a list of edges in the graph. 
            The edge notation "0->1:7" indicates that an edge connects node 0 to node 1 with weight 7.  
            You may assume a given topological order corresponding to nodes in increasing order.
     Output: The length of a longest path in the graph, followed by a longest path. 
            (If multiple longest paths exist, you may return any one.)
    """
    with open("/home/snopoff/Downloads/dataset_245_7 (1).txt", "r") as f:
        lines = f.readlines()
    source = int(lines[0].strip())
    sink = int(lines[1].strip())
    dag = defaultdict(list)
    for line in lines[2:]:
        vertex, edge_weight = line.strip().split('->')
        edge, weight = list(map(int, edge_weight.split(':')))
        dag[int(vertex)].append([edge, weight])
    length, backtrack = longest_path_in_dag(dag, source, sink)
    res = str(length) + "\n"
    v = sink
    track = "{}".format(v)
    while v != source:
        v = backtrack[v]
        track = "{}->".format(v) + track
    res += track
    with open("res.txt", "w") as f:
        f.write(res)


if __name__ == '__main__':
    forth_task()
