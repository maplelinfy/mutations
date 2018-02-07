
from multiprocessing import Pool
import pysam
import time
import re

data_path = '/data/linfengye/snpdata/'
bam_file = data_path + 'mem_pe_T20_trios_NA12878_realn_sorted.bam'
bed_file = data_path + 'hg19.bed'

base_quality_filter = 10
mapping_quality_filter = 10
read_count_filter = 1

def getBedLen(bed):
    bed_sum = 0
    for i in range(len(bed)):
        bed_sum += bed[i][2] - bed[i][1]
    return bed_sum

def readBedFile(bed_file):
    bed = []
    f = open(bed_file, 'r')
    while True:
        line = f.readline()
        if not line: break
        strr = line.split()
        bed.append([strr[0][3:], int(strr[1]), int(strr[2])])
    f.close()
    return bed

def getBedArr(bed_file, num):
    bed_list = []
    bed = readBedFile(bed_file)
    len_sum = getBedLen(bed)
    len_per = int(len_sum / num)
    bed1 = []
    need_len = len_per
    while(len(bed) != 0):
        if bed[0][2] - bed[0][1] <= need_len:
            bed1.append(bed[0])
            need_len -= bed[0][2] - bed[0][1]
            bed = bed[1:]
        else:
            bed1.append([bed[0][0], bed[0][1], bed[0][1] + need_len])
            bed[0][1] += need_len
            bed_list.append(bed1)
            bed1 = []
            need_len = len_per
        if len(bed_list) == num - 1:
            bed_list.append(bed)
            break
    for i in range(len(bed_list)):
        for j in range(len(bed_list[i])):
            bed_list[i][j][1] = str(bed_list[i][j][1])
            bed_list[i][j][2] = str(bed_list[i][j][2])
    return bed_list

def getActBase(base, base_q, cigar):
    inf = re.split(r'(\D)', cigar)
    inf = inf[0:len(inf)-1]
    act_base = ''
    act_q = []
    poi_base = 0
    poi_q = 0
    for i in range(0, len(inf), 2):
        l = int(inf[i])
        if inf[i+1] == 'M':
            act_base += base[poi_base:poi_base+l]
            act_q.extend(base_q[poi_q:poi_q+l])
            poi_base += l
            poi_q += l
        elif inf[i+1] == 'D':
            act_base += '-' * l
            act_q.extend([0]*l)
        elif inf[i+1] == 'I':
            act_base += base[poi_base:poi_base+l].lower()
            act_q.extend(base_q[poi_q:poi_q+l])
            poi_base += l
            poi_q += l
    return act_base, act_q

def recall(base, ref):
    newref = ''
    ref_poi = 0
    for i in range(len(base)):
        if base[i].isupper() or base[i] == '-':
            newref += ref[ref_poi]
            ref_poi += 1
        else:
            newref += '*'
    return newref

def updateSites(sites, chr, start_point, base, base_q, ref, cigar):
    if 'I' not in cigar and 'D' not in cigar:
        for i in range(len(base)):
            if base[i] != ref[i] and base_q[i] > base_quality_filter:
                if chr+'_'+str(start_point + i + 1)+'_'+ref[i]+'_'+base[i] in sites:
                    sites[chr+'_'+str(start_point + i + 1)+'_'+ref[i]+'_'+base[i]] += 1
                else:
                    sites[chr+'_'+str(start_point + i + 1)+'_'+ref[i]+'_'+base[i]] = 1
    else:
        base, base_q = getActBase(base, base_q, cigar)
        ref = recall(base, ref)
        pos = start_point - 1
        i = 0
        while i < len(base):
            if base[i] != '-' and ref[i] != '*':
                pos += 1
                if base[i] != ref[i]:
                    if chr + '_' + str(pos + 1) + '_' + ref[i] + '_' + base[i] in sites:
                        sites[chr + '_' + str(pos + 1) + '_' + ref[i] + '_' + base[i]] += 1
                    else:
                        sites[chr + '_' + str(pos + 1) + '_' + ref[i] + '_' + base[i]] = 1
            elif base[i] == '-':
                del_part = ref[i]
                del_len = 1
                while i + 1 < len(base) and base[i+1] == '-':
                    del_len += 1
                    i += 1
                    del_part += ref[i]
                if chr + '_' + str(pos + 1) + '_' + ref[pos-start_point] + '_-' + del_part in sites:
                    sites[chr + '_' + str(pos + 1) + '_' + ref[pos-start_point] + '_-' + del_part] += 1
                else:
                    sites[chr + '_' + str(pos + 1) + '_' + ref[pos-start_point] + '_-' + del_part] = 1
                pos += del_len
            else:
                in_part = base[i].upper()
                in_len = 1
                while i + 1 < len(base) and ref[i+1] == '*':
                    in_len += 1
                    i += 1
                    in_part += base[i].upper()
                if chr + '_' + str(pos + 1) + '_' + ref[pos-start_point] + '_+' + in_part in sites:
                    sites[chr + '_' + str(pos + 1) + '_' + ref[pos-start_point] + '_+' + in_part] += 1
                else:
                    sites[chr + '_' + str(pos + 1) + '_' + ref[pos-start_point] + '_+' + in_part] = 1
            i += 1

def getMutSites(bed_arr):
    bam = pysam.AlignmentFile(bam_file, "rb")
    sites = {}
    for i in range(len(bed_arr)):
        chr = 'chr' + str(bed_arr[i][0])
        start = int(bed_arr[i][1])
        end = int(bed_arr[i][2])
        for read in bam.fetch(chr, start-1, end):
            cigar = read.cigarstring
            read_q = int(read.mapping_quality)
            if cigar == None or read_q < mapping_quality_filter: continue
            start_point = read.pos
            base = read.query_alignment_sequence.upper()
            base_q = read.query_alignment_qualities
            ref = read.get_reference_sequence().upper()
            updateSites(sites, chr, start_point, base, base_q, ref, cigar)
    return sites

def multiprocess(bed_arr, id):
    fw = open(id + '.mutations', 'w')
    core_num = len(bed_arr)
    for i in range(core_num):
        p = Pool(core_num)
        sites_iter = p.map(getMutSites, [bed_arr[i]])
        p.close()
        p.join()
        for sites in sites_iter.keys():
            if sites_iter[sites] > read_count_filter:
                sites = sites.split('_')
                fw.write(sites[0] + '\t' + str(sites[1]) + '\t' + sites[2] + '\t' + sites[3] + '\n')
    fw.close()

def begin(core_num, id):
    bed_arr = getBedArr(bed_file, core_num)
    multiprocess(bed_arr, id)

print(time.ctime())
begin(20, 'NA12878')
print(time.ctime())
