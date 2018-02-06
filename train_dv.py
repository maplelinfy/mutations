
# -*- coding: utf-8 -*-
"""
deepvariant训练流程：
之前的encode_dv将数据都存在某个路径下snp_p,snp_n,indel_p,indel_n四个文件夹中，本程序根据数据所在路径自动获取训练以及验证数据，
由于不同类别数据比例不平衡，所有我们设定一个固定数值n，每个文件中随机选取n个数据，当n大于encode_dv中的max_num_per_file时，
即为选取所有数据进行训练。训练时我们采用chr1-18的位点，验证采用chr20-22，网络结构选取inception V3。
"""
from keras.callbacks import LearningRateScheduler
from keras.models import load_model
import numpy as np
import inceptionV3
import random
import math
import time
import os
import re

data_path = '/THL5/home/bgitj-zhuhm/mutation/data/'

def genNewY(y):
    y_new = []
    for i in range(len(y)):
        if y[i] == 0:
            y_new.append([1, 0, 0])
        elif y[i] == 1:
            y_new.append([0, 1, 0])
        else:
            y_new.append([0, 0, 1])
    return y_new

def getData(path, num):
    num0 = num1 = num2 = 0
    x_train = []
    y_train = []
    x_validate = []
    y_validate = []
    data_file = os.listdir(path)
    for i in range(len(data_file)):
        if not re.match("(.*)npy", data_file[i]): continue
        if 'y' in data_file[i].split('.')[0]: continue
        chrnum = int(data_file[i].split('_')[0][3:])
        if chrnum == 19: continue
        is_train = 1
        if chrnum < 20:
            x = x_train
            y = y_train
        else:
            x = x_validate
            y = y_validate
            is_train = 0
        y_file = path + data_file[i].split('.')[0] + '_y.npy'
        yy = np.load(y_file)
        dat = path + data_file[i]
        img = np.load(dat)
        if not is_train:
            x.extend(img)
            y.extend(yy)
            continue
        if len(img) < num:
            cap = len(dat)
        else:
            cap = num
        select_list = range(len(img))
        arr = random.sample(select_list, cap)
        for j in arr:
            x.append(img[j])
            y.append(yy[j])
            if yy[j] == 0:
                num0 += 1
            elif yy[j] == 1:
                num1 += 1
            else:
                num2 += 1
    print(num0)
    print(num1)
    print(num2)
    y_train = genNewY(y_train)
    y_validate = genNewY(y_validate)
    return x_train, y_train, x_validate, y_validate

def genTrainingData(mut_type, pos_num, neg_num):
    if mut_type == 'snp':
        pos_path = data_path + 'snp_p/'
        neg_path = data_path + 'snp_n/'
    else:
        pos_path = data_path + 'indel_p/'
        neg_path = data_path + 'indel_n/'
    x_train_pos, y_train_pos, x_validate_pos, y_validate_pos = getData(pos_path, pos_num)
    x_train_neg, y_train_neg, x_validate_neg, y_validate_neg = getData(neg_path, neg_num)
    x_train = x_train_pos + x_train_neg
    y_train = y_train_pos + y_train_neg
    x_validate = x_validate_pos + x_validate_neg
    y_validate = y_validate_pos + y_validate_neg
    x_train = np.array(x_train, dtype=np.int8)
    y_train = np.array(y_train, dtype=np.int8)
    x_validate = np.array(x_validate, dtype=np.int8)
    y_validate = np.array(y_validate, dtype=np.int8)
    np.save('x_train.npy', x_train)
    np.save('y_train.npy', y_train)
    np.save('x_validate.npy', x_validate)
    np.save('y_validate.npy', y_validate)

def validate(x_file, y_file, mut_type, cutoff):
    x_validate = np.load(x_file)
    y_validate = np.load(y_file)
    if mut_type == 'snp':
        model = load_model('snp.m')
    else:
        model = load_model('indel.m')
    y_pred = model.predict_proba(x_validate, verbose=0)
    nump = numn = tp = fp = 0
    for i in range(len(y_pred)):
        if y_pred[i][1] + y_pred[i][2] > cutoff:
            if y_validate[i][1] == 1 or y_validate[i][2] == 1:
                tp += 1
            else:
                fp += 1
        if y_validate[i][1] == 1 or y_validate[i][2] == 1:
            nump += 1
        else:
            numn += 1
    print('Positive: ' + str(nump))
    print('Negative: ' + str(numn))
    print('True positive rate: ' + str(float(tp) / nump))
    print('False positive: ' + str(fp))

def step_decay(epoch):
    initial_lr = 0.001
    epochs_drop = 2.0
    drop = 0.94
    lrate = initial_lr * math.pow(drop, math.floor((1+epoch)/epochs_drop))
    return lrate

def training(x_file, y_file, mut_type):
    x_train = np.load(x_file)
    y_train = np.load(y_file)
    model = inceptionV3.inception_v3()
    lr = LearningRateScheduler(step_decay)
    model.fit(x_train, y_train, batch_size=64, nb_epoch=30, callbacks=[lr])
    if mut_type == 'snp':
        model.save('snp.m')
    else:
        model.save('indel.m')


mut_type = 'snp'
cutoff = 0.5
pos_num_per_file = 2000
neg_num_per_file = 2000
print(time.ctime())
genTrainingData(mut_type, pos_num_per_file, neg_num_per_file)
training('x_train.npy', 'y_train.npy', mut_type)
validate('x_validate.npy', 'y_validate.npy', mut_type, cutoff)
print(time.ctime())