
ref_file = '/data/linfengye/snpdata/hg19/hg19.fasta'
data_file = '/home/linfengye/mutations/data/data_new/NA12878.mutations'
width = 221

def getPoi(data_file):
    f = open(data_file, 'r')
    poi_list = []
    while True:
        line = f.readline()
        if not line: break
        strr = line.split()
        if strr[0] == 'chrX': break
        poi_list.append([int(strr[0][3:]), int(strr[1])])
    poi_list = sorted(poi_list, key=lambda x: (x[0], x[1]))
    return poi_list

def update(ref, poi_list):
    ref_list = []
    index = 0
    chrr = poi_list[0][0]
    while chrr == poi_list[index][0]:
        start = poi_list[index][1] - (width - 1) / 2 - 1
        ref_list.append(ref[start: start+width])
        index += 1
        if index == len(poi_list): break
    if index < len(poi_list):
        poi_list = poi_list[index:]
    else:
        poi_list = []
    return ref_list, poi_list

def getRef(poi_list):
    ref = ''
    ref_list = []
    chrr = -1
    f = open(ref_file)
    while True:
        line = f.readline()
        if not line: break
        ref_now = line.split()[0].upper()
        if ref_now[0] == '>':
            if ref_now[4] == 'M': continue
            if ref_now[4] == 'X': break
            if ref != '':
                ref_list1, poi_list = update(ref, poi_list)
                ref_list.extend(ref_list1)
                if len(poi_list) == 0: break
                ref = ''
            chrr = int(ref_now[4:])
        if chrr == poi_list[0][0] and line[0] != '>':
            ref += ref_now
    if len(poi_list) != 0:
        ref_list1, poi_list = update(ref, poi_list)
        ref_list.extend(ref_list1)
    return ref_list

import time
print(time.ctime())
poi_list = getPoi(data_file)
ref_list = getRef(poi_list)
print(len(ref_list))
print(time.ctime())
