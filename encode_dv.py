#bamToPic
import numpy as np
import random
import pysam
import re

width = 221
height = 100
channels = 7

max_pixel = 254
ref_size = 5
A_value = 250
G_value = 180
T_value = 100
C_value = 30
N_value = 0
del_value = 0
insert_value = 42
base_quality_cap = 40
mapping_quality_cap = 60
positive_strand_value = 70
negative_strand_value = 240
read_support_alt_value = int(max_pixel*0.2)
read_nonsupport_alt_value = int(max_pixel*1)
base_match_value = int(max_pixel*1)
base_mismatch_value = int(max_pixel*0.6)

base_quality_filter = 10
mapping_quality_filter = 10
read_count_filter = 1
max_num_per_file = 8000

def randomIndex(n, hei):
    select_list = range(n)
    arr = random.sample(select_list, hei)
    arr.sort()
    return arr

def getActBase(base, base_q, inf):
    inf = re.split(r'(\D)', inf)
    inf = inf[0:len(inf)-1]
    act_base = []
    act_q = []
    base = list(base)
    poi_base = 0
    poi_q = 0
    for i in range(0, len(inf), 2):
        l = int(inf[i])
        if inf[i+1] == 'M':
            act_base.extend(base[poi_base:poi_base+l])
            act_q.extend(base_q[poi_q:poi_q+l])
            poi_base += l
            poi_q += l
        elif inf[i+1] == 'D':
            act_base.extend(['-'] * l)
            act_q.extend([0]*l)
        elif inf[i+1] == 'I':
            act_base.extend(['*'] * l)
            act_q.extend(base_q[poi_q:poi_q+l])
            poi_base += l
            poi_q += l
    act_base = ''.join(act_base)
    return act_base, act_q

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
    elif base == '*':
        return insert_value
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
                            for p in range(channels):
                                arr[k][t][p] = arr[k][t-1][p]
                        for p in range(channels):
                            arr[k][i][p] = 0
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

def ifMut(base, ref, k, start_point, cigar, alt):
    inf = re.split(r'(\D)', cigar)
    inf = inf[0:len(inf)-1]
    if 'I' in cigar or 'D' in cigar:
        diff = k - start_point + 1
        summ = 0
        for i in range(0, len(inf), 2):
            l = int(inf[i])
            if inf[i+1] == 'M' or inf[i+1] == 'D':
                summ += l
                if summ == diff:
                    return 1
        return 0
    else:
        if base[k-start_point] != ref[k-start_point] and base[k-start_point] == alt:
            return 1
        else:
            return 0

def updateRef(arr_ref, ref, start_point, pic_start, pic_end):
    if start_point < pic_start:
        arr_start = 0
        ref_start = pic_start - start_point
    else:
        arr_start = start_point - pic_start
        ref_start = 0
    if start_point + len(ref) - 1 >= pic_end:
        arr_end = len(arr_ref)
    else:
        arr_end = arr_start + len(ref)
    for i in range(arr_start, arr_end):
        arr_ref[i] = getBaseValue(ref[ref_start])
        ref_start += 1
        if ref_start >= len(ref):
            break
    return arr_ref

def buildRef(ref):
    refnew = []
    for i in range(width):
        if ref[i] != 0:
            refnew.append([ref[i], max_pixel, max_pixel, positive_strand_value, max_pixel, max_pixel, 0])
        else:
            refnew.append([0]*channels)
    ref = [refnew] * ref_size
    return ref

def addCigar(read_arr,  start_point, pic_start, pic_end, cigar):
    cigar_arr = [0] * width
    cigar_list = []
    inf = re.split(r'(\D)', cigar)
    inf = inf[0:len(inf) - 1]
    for i in range(0, len(inf), 2):
        if inf[i+1] == 'M' or inf[i+1] == 'I':
            cigar_list.extend([int(inf[i])] * int(inf[i]))
        else:
            cigar_list.extend([0] * int(inf[i]))
    if start_point < pic_start:
        arr_start = 0
        cigar_start = pic_start - start_point
    else:
        arr_start = start_point - pic_start
        cigar_start = 0
    if start_point + len(cigar_list) - 1 >= pic_end:
        arr_end = len(cigar_arr)
    else:
        arr_end = arr_start + len(cigar_list)
    for i in range(arr_start, arr_end):
        cigar_arr[i] = cigar_list[cigar_start]
        cigar_start += 1
        if cigar_start >= len(cigar_list):
            break
    for i in range(len(read_arr)):
        read_arr[i][6] = cigar_arr[i]
    return read_arr

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
                for t in range(channels):
                    arr[i][j][t] = 0
        else:
            for j in range(enp-toz, enp+1):
                for t in range(channels):
                    arr[i][j][t] = 0
    return arr

def savefile(path, list_name, y, data, id, chr):
    data = np.array(data, dtype=np.int8)
    np.save(path + chr + '_' + str(id) + '.npy', data)
    y = np.array(y, dtype=np.int8)
    np.save(path + chr + '_' + str(id) + '_y.npy', y)
    file_name = path + chr + '_' + str(id) + '.list'
    wf = open(file_name, 'w')
    for i in range(len(list_name)):
        wf.write(list_name[i] + '\n')
    wf.close()

def genReadArr(pic_start, pic_end, start_point, k, base, base_q, read_q, ref, strand, cigar, alt):
    read_arr = []
    for i in range(pic_start, pic_end):
        poi = [0] * channels
        t = i - start_point
        if t < 0 or t > len(base) - 1:
            read_arr.append(poi)
            continue
        poi[0] = getBaseValue(base[t])
        poi[1] = int(min(int(base_q[t]), base_quality_cap) * max_pixel / base_quality_cap)
        if base[t] != '-':
            poi[2] = int(min(read_q, mapping_quality_cap) * max_pixel / mapping_quality_cap)
        if strand:
            poi[3] = negative_strand_value
        else:
            poi[3] = positive_strand_value
        if ifMut(base, ref, k, start_point, cigar, alt):
            poi[4] = read_nonsupport_alt_value
        else:
            poi[4] = read_support_alt_value
        if base[t] != ref[t] or base[t] == '*' or base[t] == '-':
            poi[5] = base_match_value
        else:
            poi[5] = base_mismatch_value
        read_arr.append(poi)
    read_arr = addCigar(read_arr, start_point, pic_start, pic_end, cigar)
    return read_arr

def genReadMatrix(bamfile, chr, k, alt):
    arr = arrtemp = []
    arr_ref = [0] * width
    for read in bamfile.fetch(chr, k, k + 1):
        start_point = read.pos
        base_q = read.query_alignment_qualities
        if base_q[k-start_point] < base_quality_filter: continue
        read_q = int(read.mapping_quality)
        if read.cigarstring == None or read_q < mapping_quality_filter: continue
        base = read.query_alignment_sequence
        base = base.upper()
        ref = read.get_reference_sequence()
        ref = ref.upper()
        pic_start = k - int(width/2)
        if width % 2 == 0:
            pic_end = k + int(width/2)
        else:
            pic_end = k + int(width/2) + 1
        arr_ref = updateRef(arr_ref, ref, start_point, pic_start, pic_end)
        strand = read.is_reverse
        cigar = read.cigarstring
        if 'I' in cigar or 'D' in cigar:
            base, base_q = getActBase(base, base_q, cigar)
        if 'I' in cigar:
            ref = recall(base, ref)
        read_arr = genReadArr(pic_start, pic_end, start_point, k, base, base_q, read_q, ref, strand, cigar, alt)
        arrtemp.append(read_arr)
    arr_ref = buildRef(arr_ref)
    depth = len(arrtemp)
    height_reads = height - ref_size
    if depth > height_reads:
        index = randomIndex(depth, height_reads)
        for i in range(height_reads):
            arr.append(arrtemp[index[i]])
        arr = [arr_ref, arr]
        arr = adjust(arr)
    elif depth == height_reads:
        arr = [arr_ref, arrtemp]
        arr = adjust(arr)
    elif depth < height_reads and depth >= read_count_filter:
        arr = [arr_ref, arrtemp]
        arr = adjust(arr)
        for i in range(height_reads - depth):
            arr.append([[0] * channels] * width)
    return arr

def code(bam_file, file, path):
    id = 1
    chr_pre = 'chr'
    name = data = y = []
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    sf = open(file, 'r')
    while True:
        line = sf.readline()
        if not line: break
        if line[0] == '#': continue
        strr = line.split()
        chr = strr[0]
        if chr_pre == 'chr':
            chr_pre = chr
        if chr != chr_pre:
            savefile(path, name, y, data, id, chr_pre)
            id = 1
            name = data = y = []
        k = int(strr[1]) - 1
        ref = strr[2]
        alt = strr[3]
        tp = int(strr[-1])
        arr = genReadMatrix(bamfile, chr, k, alt)
        if len(arr) < 1:
            continue
        data.append(arr)
        y.append(tp)
        sample_id = chr + ':' + str(k+1) + ':' + ref + ':' + alt
        name.append(sample_id)
        chr_pre = chr
        if len(data) == max_num_per_file:
            savefile(path, name, y, data, id, chr_pre)
            id += 1
            name = data = y = []
    savefile(path, name, y, data, id, chr_pre)
    sf.close()


data_path = '/THL4/home/bgi_linfengye/mutations/data/'
bam_file = data_path + 'mem_pe_T20_trios_NA12878_realn_sorted.bam'
snp_positive_path = data_path + 'snp_p/'
snp_negative_path = data_path + 'snp_n/'
indel_positive_path = data_path + 'indel_p/'
indel_negative_path = data_path + 'indel_n/'
snp_positive = data_path + 'bed.na878.mut.snp.nofilt.TP'
snp_negative = data_path + 'bed_na878.mut.snp.nofilt.FP'
indel_positive = data_path + 'bed.na878.mut.indel.nofilt.TP'
indel_negative = data_path + 'bed.na878.mut.indel.nofilt.FP'

import time
print(time.ctime())
code(bam_file, snp_positive, snp_positive_path)
print(time.ctime())
