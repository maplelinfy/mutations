
# -*- coding: utf-8 -*-
"""
参考链接 https://blog.csdn.net/g15827636417/article/details/52818815

计算稀疏矩阵mat和mat转置的乘积：
get_non_zero_arr函数是根据矩阵mat获取如下四个信息
nonzero_arr为mat中所有非零元素的list，每个元素依次为横坐标，纵坐标，值
nonzero_arr_tran为mat转置的非零元素list
nonzero_num_col为mat转置的每一行中非零元素的个数
nonzero_index_col为mat转置的每一行中第一个非零元素在nonzero_arr_tran中的位置

设矩阵规格为m*n,非零元素个数为p
获取信息的时间复杂度为O(m*n)，乘法中对结果矩阵初始化时间复杂度为O(m*m)，做乘法过程中的时间复杂度为O(p*p/n)
"""

import numpy as np

def get_non_zero_arr(mat):
    row = mat.shape[0]
    col = mat.shape[1]
    nonzero_arr = []
    nonzero_arr_tran = []
    nonzero_num_col = np.zeros(col, dtype=np.int32)
    nonzero_index_col = np.zeros(col, dtype=np.int32)
    for i in range(row):
        for j in range(col):
            if mat[i][j] != 0:
                nonzero_arr.append([i, j, mat[i][j]])
    for j in range(col):
        for i in range(row):
            if mat[i][j] != 0:
                nonzero_arr_tran.append([j, i, mat[i][j]])
                nonzero_num_col[j] += 1
        if j != col - 1:
            nonzero_index_col[j+1] = nonzero_index_col[j] + nonzero_num_col[j]
    return nonzero_arr, nonzero_arr_tran, nonzero_num_col, nonzero_index_col

def multiplication_of_sparse_matrix(mat):
    row = mat.shape[0]
    col = mat.shape[1]
    nonzero_arr, nonzero_arr_tran , nonzero_num_col, nonzero_index_col = get_non_zero_arr(mat)
    nonzero_num = len(nonzero_arr_tran)
    res = np.zeros(shape=(row, row))
    p = 0
    for i in range(row):
        while(p < nonzero_num and i == nonzero_arr[p][0]):
            k = nonzero_arr[p][1]
            if k < col - 1:
                t = nonzero_index_col[k+1]
            else:
                t = nonzero_num
            for q in range(nonzero_index_col[k], t):
                j = nonzero_arr_tran[q][1]
                res[i][j] += nonzero_arr[p][2] * nonzero_arr_tran[q][2]
            p += 1
        print(res[i])
    return res

# test
mat = np.ones(shape=(3, 2))
res = multiplication_of_sparse_matrix(mat)
print(res)