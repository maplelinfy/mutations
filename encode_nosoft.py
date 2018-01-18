#################################
# snp and indel without softclip
#################################
import numpy as np
import random
import pysam
import re

width = 200
height = 60
max_num = 20000
bp_len = 90

A_value = 50
C_value = 55
G_value = 60
T_value = 65
N_value = 40
del_value = 1
insert_value = 80
match_value = 10
quality_filter = 10
read_count_filter = 1

def randomIndex(n, hei):
    select_list = range(n)
    arr = random.sample(select_list, hei)
    arr.sort()
    return arr

def getActBase(base, read_q, inf):
    inf = re.split(r'(\D)', inf)
    inf = inf[0:len(inf)-1]
    actbase = []
    act_q = []
    base = list(base)
    poi_base = 0
    poi_q = 0
    for i in range(0,len(inf),2):
        l = int(inf[i])
        if inf[i+1] == 'M':
            actbase.extend(base[poi_base:poi_base+l])
            act_q.extend(read_q[poi_q:poi_q+l])
            poi_base += l
            poi_q += l
        elif inf[i+1] == 'D':
            actbase.extend(['-']*l)
            act_q.extend([0]*l)
        elif inf[i+1] == 'I':
            actbase.extend(['*'] * l)
            act_q.extend(read_q[poi_q:poi_q + l])
            poi_base += l
            poi_q += l
    actbase = ''.join(actbase)
    return actbase, act_q

def getBaseValue(base):
    if base == 'A':
        return A_value
    elif base == 'C':
        return C_value
    elif base == 'G':
        return G_value
    elif base == 'T':
        return T_value
    elif base == '-':
        return del_value
    else:
        return N_value

def adjust(arr):
    for i in range(width):
        for j in range(len(arr)):
            if arr[j][i][0] == insert_value:
                for k in range(len(arr)):
                    if arr[k][i][0] != insert_value:
                        a = list(range(i+1, width))
                        a.reverse()
                        for t in a:
                            arr[k][t][0] = arr[k][t-1][0]
                            arr[k][t][1] = arr[k][t-1][1]
                            arr[k][t][2] = arr[k][t-1][2]
                        arr[k][i][0] = 0
                        arr[k][i][1] = 0
                        arr[k][i][2] = 0
                break
    return arr

def recall(base, ref):
    newref = ''
    ref_poi = 0
    for i in range(len(base)):
        if base[i] != '*':
            newref += ref[ref_poi]
            ref_poi += 1
        else:
            newref += '*'
    return newref

def intercept(arr, length):
    for i in range(len(arr)):
        stp = 0
        enp = width - 1
        while arr[i][stp][0] == 0:
            stp += 1
            if stp == width: break
        while arr[i][enp][0] == 0:
            enp -= 1
            if enp < 0: break
        read_len = enp - stp + 1
        if read_len <= length: continue
        k = int(width/2) - 1
        toz = read_len - length
        if k - stp > enp - k:
            for j in range(stp, stp+toz+1):
                arr[i][j][0] = 0
                arr[i][j][1] = 0
                arr[i][j][2] = 0
        else:
            for j in range(enp-toz, enp+1):
                arr[i][j][0] = 0
                arr[i][j][1] = 0
                arr[i][j][2] = 0
    return arr

def savefile(list_name, file_name):
    wf = open(file_name, 'w')
    for i in range(len(list_name)):
        wf.write(list_name[i] + '\n')
    wf.close()

def if_mut(base, ref, k, start_point, cigar):
    inf = re.split(r'(\D)', cigar)
    inf = inf[0:len(inf)-1]
    diff = k - start_point + 1
    summ = 0
    if 'I' in cigar or 'D' in cigar:
        for i in range(0, len(inf), 2):
            l = int(inf[i])
            if inf[i+1] == 'M' or inf[i+1] == 'D':
                summ += l
                if summ == diff:
                    return 1
        return 0
    else:
        if base[k - start_point] != ref[k - start_point]:
            return 1
        else:
            return 0

def code(bam_file, file, path):
    sf = open(file, 'r')
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    data = []
    y = []
    name = []
    id = 1
    count = 0
    chr_pre = 'chr'
    while True:
        line = sf.readline()
        if not line: break
        if line[0] == '#': continue
        stri = line.split()
        chr = stri[0]
        if count == 0:
            chr_pre = chr
        if chr != chr_pre:
            data = np.array(data, dtype=np.int8)
            np.save(path+chr_pre+'_'+str(id)+'.npy', data)
            y = np.array(y, dtype=np.int8)
            np.save(path+chr_pre+'_'+str(id)+'_y.npy', y)
            y = []
            savefile(name, path+chr_pre+'_'+str(id)+'.list')
            data = []
            name = []
            id = 1
        k = int(stri[1]) - 1
        tp = int(stri[-1])
        arr = []
        arrtemp_p = []
        arrtemp_n = []
        for read in bamfile.fetch(chr, k, k+1):
            start_point = read.pos
            read_q = int(read.mapping_quality)
            if read.cigarstring == None or read_q < quality_filter: continue
            base = read.query_alignment_sequence
            base = base.upper()
            base_q = read.query_alignment_qualities
            ref = read.get_reference_sequence()
            ref = ref.upper()
            if 'I' in read.cigarstring or 'D' in read.cigarstring:
                base, base_q = getActBase(base, base_q, read.cigarstring)
                if start_point + len(base) <= k: continue
            if 'I' in read.cigarstring:
                ref = recall(base, ref)
            arrtemp = arrtemp_n
            if if_mut(base, ref, k, start_point, read.cigarstring):
                arrtemp = arrtemp_p
            readarr = []
            for i in range(k-int(width/2)+1, k+int(width/2)+1):
                poi = [0, 0, 0]
                t = i - start_point
                if t < 0 or t > len(base)-1:
                    readarr.append(poi)
                    continue
                if base[t] != ref[t]:
                    poi[0] = getBaseValue(base[t])
                elif base[t] == '*':
                    poi[0] = insert_value
                else:
                    poi[0] = match_value
                poi[1] = int(base_q[t])
                if base[t] != '-':
                    poi[2] = read_q
                readarr.append(poi)
            arrtemp.append(readarr)
        arrtemp = []
        if len(arrtemp_p) != 0:
            for i in range(len(arrtemp_p)):
                arrtemp.append(arrtemp_p[i])
        if len(arrtemp_n) != 0:
            for i in range(len(arrtemp_n)):
                arrtemp.append(arrtemp_n[i])
        num = len(arrtemp)
        if num > height:
            index = randomIndex(num, height)
            for i in range(height):
                arr.append(arrtemp[index[i]])
            arr = adjust(arr)
        elif num == height:
            arr = adjust(arrtemp)
        elif num < height and num >= read_count_filter:
            arr = adjust(arrtemp)
            for i in range(height - num):
                arr.append([[0]*3]*width)
        if len(arr) < 1:
            continue
        count += 1
        chr_pre = chr
        data.append(arr)
        y.append(tp)
        sampleID = chr + '_' + str(k+1)
        name.append(sampleID)
        if len(data) == max_num:
            data = np.array(data, dtype=np.int8)
            np.save(path+chr+'_'+str(id)+'.npy', data)
            y = np.array(y, dtype=np.int8)
            np.save(path + chr_pre + '_' + str(id) + '_y.npy', y)
            y = []
            savefile(name, path+chr+'_'+str(id)+'.list')
            id += 1
            data = []
            name = []
    data = np.array(data, dtype=np.int8)
    np.save(path+chr_pre+'_'+str(id)+'.npy', data)
    y = np.array(y, dtype=np.int8)
    np.save(path + chr_pre + '_' + str(id) + '_y.npy', y)
    savefile(name, path+chr_pre+'_'+str(id)+'.list')
    sf.close()


data_path = '/THL4/home/bgi_linfengye/mutations/data/'
snp_positive_path = data_path + 'snp_p/'
snp_negative_path = data_path + 'snp_n/'
indel_positive_path = data_path + 'indel_p/'
indel_negative_path = data_path + 'indel_n/'
snp_positive = data_path + 'snp_pos_filt.list'
snp_negative = data_path + 'bed_na878.mut.snp.filt.FP'
indel_positive = data_path + 'indel_pos.list'
indel_negative = data_path + 'bed_na878.mut.indel.nofilt.FP'
bam_file = data_path + 'mem_pe_T20_trios_NA12878_realn_sorted.bam'

import time
print(time.ctime())
code(bam_file, snp_positive, snp_positive_path)
print(time.ctime())