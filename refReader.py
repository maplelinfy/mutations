# Read the reference by fixed lines
ref_file = '/data/linfengye/snpdata/hg19/hg19.fasta'
data_file = '/home/linfengye/mutations/data/data_new/NA12878.mutations'
width = 221
line_num = 1000000

def update(ref, poi_list, chrr, ref_start, ref_end):
    flag = 0
    size_list = []
    ref_list1 = []
    poi_index = 0
    while poi_list[poi_index][0] == chrr and \
            poi_list[poi_index][1] - 1 - (width - 1) / 2 >= ref_start and \
            poi_list[poi_index][1] - 1 - (width - 1) / 2 <= ref_end:
        index1 = poi_list[poi_index][1] - 1 - (width - 1) / 2 - ref_start
        if poi_list[poi_index][1] - 1 + (width - 1) / 2 <= ref_end:
            ref_list1.append(ref[index1: index1+width])
        else:
            ref_list1.append(ref[index1:])
            size_list.append(poi_list[poi_index][1] - 1 + (width - 1) / 2 - ref_end)
            flag = 1
        poi_index += 1
        if poi_index == len(poi_list): break
    if poi_index < len(poi_list):
        poi_list = poi_list[poi_index:]
    else:
        poi_list = []
    return flag, size_list, poi_list, ref_list1

def getRef(poi_list):
    ref = ''
    chrr = -1
    ref_start = 0
    ref_end = -1
    line_now = 0
    ref_list = []
    flag = 0
    size_list = []

    f = open(ref_file, 'r')
    while True:
        line = f.readline()
        if not line: break
        ref_now = line.split()[0].upper()
        if ref_now[0] == '>':
            if ref_now[4] == 'X': break
            if ref_now[4] == 'M': continue
            if ref != '':
                flag, size_list, poi_list, ref_list1 = update(ref, poi_list, chrr, ref_start, ref_end)
                ref_list.extend(ref_list1)
                if len(poi_list) == 0: break
                ref_start = 0
                ref_end = -1
                line_now = 0
                ref = ''
            chrr = int(ref_now[4:])
        if chrr == poi_list[0][0] and line[0] != '>':
            if flag:
                for i in range(len(size_list)):
                    if size_list[i] <= len(ref_now):
                        ref_list[i - len(size_list)] += ref_now[0: size_list[i]]
                        size_list[i] = 0
                    else:
                        ref_list[i-len(size_list)] += ref_now
                        size_list[i] -= len(ref_now)
                if sum(size_list) == 0:
                    flag = 0
                    if len(poi_list) == 0: break
            if line_now < line_num:
                line_now += 1
            else:
                flag, size_list, poi_list, ref_list1 = update(ref, poi_list, chrr, ref_start, ref_end)
                ref_list.extend(ref_list1)
                if len(poi_list) == 0 and flag == 0: break
                line_now = 1
                ref_start = ref_end + 1
                ref = ''
            ref += ref_now
            ref_end += len(ref_now)
    if len(poi_list) != 0:
        flag, size_list, poi_list, ref_list1 = update(ref, poi_list, chrr, ref_start, ref_end)
        ref_list.extend(ref_list1)
    return ref_list

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

import time
print(time.ctime())
poi_list = getPoi(data_file)
ref_list = getRef(poi_list)
print(len(ref_list))
print(time.ctime())
