
# -*- coding: utf-8 -*-
"""
为候选突变位点添加ref流程：
输入：候选突变位点文件，ref文件，输出文件名
输出：添加ref后的候选突变位点文件
"""

ref_file = '/data/linfengye/snpdata/hg19/hg19.fasta'
data_file = '/home/linfengye/mutations/dv/NA12878.mutations'
out_file = '/home/linfengye/mutations/dv/NA12878.mutations2'
width = 221

def getPoi(data_file):
    f = open(data_file, 'r')
    poi_list = []
    while True:
        line = f.readline()
        if not line: break
        strr = line.split()
        poi_list.append([int(strr[0][3:]), int(strr[1]), strr[2], strr[3]])
    poi_list = sorted(poi_list, key=lambda x: (x[0], x[1]))
    return poi_list

def addDelRef(pos, alt, ref):
    del_len = len(alt) - 1
    del_ref = ref[pos: pos+del_len]
    return del_ref

def update(ref, poi_list, f):
    index = 0
    chrr = poi_list[0][0]
    while chrr == poi_list[index][0]:
        if len(poi_list[index][3]) > 1:
            if poi_list[index][3][1] == '-':
                poi_list[index][3] = addDelRef(poi_list[index][1], poi_list[index][3], ref)
        start = poi_list[index][1] - (width - 1) / 2 - 1
        f.write('chr' + str(poi_list[index][0]) + '\t' + str(poi_list[index][1]) + '\t' +
                poi_list[index][2] + '\t' + poi_list[index][3] + '\t' + ref[start: start+width] + '\n')
        index += 1
        if index == len(poi_list): break
    if index < len(poi_list):
        poi_list = poi_list[index:]
    else:
        poi_list = []
    return poi_list

def getRef(poi_list):
    ref = ''
    chrr = -1
    f = open(ref_file, 'r')
    fw = open(out_file, 'w')
    while True:
        line = f.readline()
        if not line: break
        ref_now = line.split()[0].upper()
        if ref_now[0] == '>':
            if ref_now[4] == 'M' or ref_now[4] == 'X': continue
            if ref != '':
                poi_list = update(ref, poi_list, fw)
                if len(poi_list) == 0: break
                ref = ''
            chrr = int(ref_now[4:])
        if chrr == poi_list[0][0] and line[0] != '>':
            ref += ref_now
    if len(poi_list) != 0:
        update(ref, poi_list, fw)
    f.close()
    fw.close()

import time
print(time.ctime())
poi_list = getPoi(data_file)
getRef(poi_list)
print(time.ctime())
