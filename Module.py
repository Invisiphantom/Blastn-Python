import numpy as np


def cal_score(str1, str2, matchScore=2, mismatchScore=-3, gapScore=-4):
    """计算比对分数矩阵 (动态规划)"""
    seq1, seq2 = list("-" + str1), list("-" + str2)
    n1, n2 = len(seq1), len(seq2)
    matrix = np.zeros((n1, n2), dtype=int)

    # 全局比对, 初始化首行列为留空值
    for i in range(1, n1):
        matrix[i][0] = gapScore * i
    for j in range(1, n2):
        matrix[0][j] = gapScore * j

    # 动态规划搜索全部路径
    for i in range(1, n1):
        for j in range(1, n2):
            lScore = matrix[i][j - 1] + gapScore  # 左侧seq1留空
            uScore = matrix[i - 1][j] + gapScore  # 上侧seq2留空
            # 对角线尝试匹配
            if seq1[i] == seq2[j]:
                diagScore = matrix[i - 1][j - 1] + matchScore
            else:
                diagScore = matrix[i - 1][j - 1] + mismatchScore
            matrix[i][j] = max(lScore, uScore, diagScore)

    return matrix


def trace_back(str1, str2, mat, matchScore=2, mismatchScore=-3, gapScore=-4):
    """回溯记录所有路径"""
    seq1, seq2 = list("-" + str1), list("-" + str2)
    n1, n2 = len(seq1), len(seq2)

    # 延申比对, 从末尾行列的最大值开始回溯
    maxI, maxJ = n1 - 1, n2 - 1
    maxVal = mat[maxI][maxJ]

    # 遍历最后一行
    for j in range(n2 - 1, 0, -1):
        if mat[n1 - 1][j] > maxVal:
            maxI, maxJ = n1 - 1, j
            maxVal = mat[n1 - 1][j]

    # 遍历最后一列
    for i in range(n1 - 1, 0, -1):
        if mat[i][n2 - 1] > maxVal:
            maxI, maxJ = i, n2 - 1
            maxVal = mat[i][n2 - 1]

    # 设置路径终点
    mStack = [(maxI, maxJ)]

    path = []
    nborStack = []
    while True:
        # 全局比对, 搜索直到起点结束
        if mStack[-1] == (0, 0):
            path = mStack.copy()
            mStack.pop()
            break

        # 记录当前位置 有效的路径方向
        if len(mStack) == len(nborStack) + 1:
            nborStack.append([])
            row = mStack[-1][0]
            col = mStack[-1][1]

            # 来自左侧
            if mat[row][col] == mat[row][col - 1] + gapScore:
                nborStack[-1].append((row, col - 1))

            # 来自上侧
            if mat[row][col] == mat[row - 1][col] + gapScore:
                nborStack[-1].append((row - 1, col))

            # 来自对角线
            if seq1[row] == seq2[col]:
                if mat[row][col] == mat[row - 1][col - 1] + matchScore:
                    nborStack[-1].append((row - 1, col - 1))
            else:
                if mat[row][col] == mat[row - 1][col - 1] + mismatchScore:
                    nborStack[-1].append((row - 1, col - 1))

        # 还存在有效方向, 继续向下走
        elif nborStack[-1] != []:
            mStack.append(nborStack[-1].pop())

        # 已经无有效方向, 则回溯
        elif nborStack[-1] == []:
            nborStack.pop()
            mStack.pop()

    return maxI, maxJ, maxVal, path


def format_path(str1, str2, path, is_forward):
    """打印路径对应的配对方式"""
    seq1, seq2 = list("-" + str1), list("-" + str2)

    printSeq1 = []
    printSeq2 = []
    printConnect = []
    for i in range(1, len(path)):
        # 来自左侧, seq1留空
        if path[i][0] == path[i - 1][0]:
            printSeq1.append("-")
            printSeq2.append(seq2[path[i - 1][1]])
            printConnect.append(" ")
        # 来自上侧, seq2留空
        elif path[i][1] == path[i - 1][1]:
            printSeq1.append(seq1[path[i - 1][0]])
            printSeq2.append("-")
            printConnect.append(" ")
        # 来自对角线
        else:
            printSeq1.append(seq1[path[i - 1][0]])
            printSeq2.append(seq2[path[i - 1][1]])
            if seq1[path[i - 1][0]] == seq2[path[i - 1][1]]:
                printConnect.append("|")
            else:
                printConnect.append("*")

    if is_forward:
        printSeq1.reverse()
        printSeq2.reverse()
        printConnect.reverse()
    return "".join(printSeq1) + "\n" + "".join(printConnect) + "\n" + "".join(printSeq2) + "\n\n"


