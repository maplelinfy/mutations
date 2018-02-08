
gold_snp_file = '/data/linfengye/snpdata/NA12878.wgs.illumina_platinum.20140404.snps_v2.vcf-filt'
gold_indel_file = '/data/linfengye/snpdata/NA12878.wgs.illumina_platinum.20140404.indels_v2.vcf-TP'
bed_file = '/data/linfengye/snpdata/NA12878/v3.3.2/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed'
datapath = '/home/linfengye/mutations/dv/'
mut_file = datapath + 'NA12878.mutations2'
snp_pos_file = datapath + 'na878.mut.snp.TP'
snp_neg_file = datapath + 'na878.mut.snp.FP'
indel_pos_file = datapath + 'na878.mut.indel.TP'
indel_neg_file = datapath + 'na878.mut.indel.FP'
bed_snp_pos_file = datapath+'bed.na878.mut.snp.TP'
bed_snp_neg_file = datapath + 'bed.na878.mut.snp.FP'
bed_indel_pos_file = datapath + 'bed.na878.mut.indel.TP'
bed_indel_neg_file = datapath + 'bed.na878.mut.indel.FP'


def splitSNPandINDEL(mut_file):
    f = open(mut_file, 'r')
    snp = []
    indel = []
    while True:
        line = f.readline()
        if not line: break
        strr = line.split()
        if strr[0] == 'chrX' or strr[0] == 'chrY' or strr[0] == 'chrM': continue
        temp = snp
        if len(strr[2]) != len(strr[3]):
            temp = indel
        strr[0] = int(strr[0][3:])
        strr[1] = int(strr[1])
        temp.append(strr)
    f.close()
    snp = sorted(snp, key=lambda x: (x[0], x[1]))
    indel = sorted(indel, key=lambda x: (x[0], x[1]))
    return snp, indel

def getBed(bed_file):
    arr_bed = []
    f = open(bed_file, 'r')
    while True:
        line = f.readline()
        if not line: break
        strr = line.split()
        if strr[0] == 'X': break
        arr_bed.append([int(strr[0]), int(strr[1]), int(strr[2])])
    f.close()
    return arr_bed

def getGoldArr(gold_file):
    arr = []
    f = open(gold_file, 'r')
    while True:
        line = f.readline()
        if not line: break
        strr = line.split()
        if strr[0] == 'chrX': break
        if strr[-1][0] == strr[-1][2]:
            label = '2'
        else:
            label = '1'
        arr.append([int(strr[0][3:]), int(strr[1]), label])
    f.close()
    return arr

def candiSplit(arr_a, gold_file, out_pos_file, out_neg_file):
    arr_b = getGoldArr(gold_file)
    f_pos = open(out_pos_file, 'w')
    f_neg = open(out_neg_file, 'w')
    index = 0
    for i in range(len(arr_a)):
        chrr = arr_a[i][0]
        pos = arr_a[i][1]
        while (chrr == arr_b[index][0] and pos > arr_b[index][1]) or (chrr > arr_b[index][0]):
            index += 1
            if index == len(arr_b):
                break
        if index < len(arr_b):
            if chrr != arr_b[index][0] or pos != arr_b[index][1]:
                f_neg.write('chr' + str(arr_a[i][0]) + '\t' + str(arr_a[i][1]) + '\t')
                for j in range(2, len(arr_a[i])):
                    f_neg.write(arr_a[i][j] + '\t')
                f_neg.write('0' + '\n')
                continue
            else:
                f_pos.write('chr' + str(arr_a[i][0]) + '\t' + str(arr_a[i][1]) + '\t')
                for j in range(2, len(arr_a[i])):
                    f_pos.write(arr_a[i][j] + '\t')
                f_pos.write(arr_b[index][2] + '\n')
        else:
            break

def filtbed(file, arr_b, outfile):
    f = open(file, 'r')
    arr_a = []
    while True:
        line = f.readline()
        if not line: break
        strr = line.split()
        strr[0] = int(strr[0][3:])
        strr[1] = int(strr[1])
        arr_a.append(strr)
    fw = open(outfile, 'w')
    index = 0
    for i in range(len(arr_a)):
        chrr = arr_a[i][0]
        pos = arr_a[i][1]
        while (chrr == arr_b[index][0] and pos > arr_b[index][2]) or (chrr > arr_b[index][0]):
            index += 1
            if index == len(arr_b):
                break
        if index < len(arr_b):
            if chrr != arr_b[index][0] or pos < arr_b[index][1]:
                continue
            else:
                fw.write('chr' + str(arr_a[i][0]) + '\t' + str(arr_a[i][1]) + '\t')
                for j in range(2, len(arr_a[i])):
                    fw.write(arr_a[i][j] + '\t')
                fw.write('\n')
        else:
            break


def begin():
    snp, indel = splitSNPandINDEL(mut_file)
    candiSplit(snp, gold_snp_file, snp_pos_file, snp_neg_file)
    candiSplit(indel, gold_indel_file, indel_pos_file, indel_neg_file)
    arr_bed = getBed(bed_file)
    filtbed(snp_pos_file, arr_bed, bed_snp_pos_file)
    filtbed(snp_neg_file, arr_bed, bed_snp_neg_file)
    filtbed(indel_pos_file, arr_bed, bed_indel_pos_file)
    filtbed(indel_neg_file, arr_bed, bed_indel_neg_file)

begin()