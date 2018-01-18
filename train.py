from keras.layers.convolutional import Conv2D, MaxPooling2D
from keras.layers.normalization import BatchNormalization
from keras.layers.core import Dense, Dropout, Flatten
from sklearn.model_selection import StratifiedKFold
from keras.callbacks import ReduceLROnPlateau
from keras.callbacks import ModelCheckpoint
from keras.callbacks import EarlyStopping
from keras.models import Sequential
from keras.models import load_model
from keras.optimizers import SGD
import numpy as np
import time
import os
import re

img_rows = 60
img_cols = 200
data_path = '/THL5/home/bgitj-zhuhm/mutation/data/'

def getname(data_file):
    file_name = data_file.split('.')[0] + '.list'
    name = []
    f = open(file_name, 'r')
    while True:
        line = f.readline()
        if not line: break
        name.append(line.split()[0])
    return name

def validate(type):
    if type == 'indel':
        positive_path = data_path + 'indel_p/'
        negative_path = data_path + 'indel_n/'
    else:
        positive_path = data_path + 'snp_p/'
        negative_path = data_path + 'snp_n/'
    sump = 0; sumn = 0; tp = 0; fp = 0
    f = open(type+'.xls', 'w')
    pos = os.listdir(positive_path)
    for i in range(len(pos)):
        if not re.match("(.*)npy", pos[i]): continue
        if 'y' in pos[i].split('.')[0]: continue
        chrnum = int(pos[i].split('_')[0][3:])
        if chrnum % 2 == 0 or chrnum == 7 or chrnum == 19: continue
        dat = positive_path + pos[i]
        name = getname(dat)
        img = np.load(dat)
        if len(img) > 0:
            if type == 'indel':
                model = load_model('indel.m')
            else:
                model = load_model('snp.m')
            # y_pred = model.predict_classes(img, verbose=0).T[0]
            # sump += len(y_pred)
            # for j in range(len(y_pred)):
            #     if y_pred[j] == 1:
            #         tp += 1
            #     else:
            #         f.write(name[j] + '\t' + str(1) + '\n')
            y_pred = model.predict_proba(img, verbose=0)
            sump += len(y_pred)
            for j in range(len(y_pred)):
                if y_pred[j][1] + y_pred[j][2] > 0.5:
                    tp += 1
                else:
                    f.write(name[j] + '\t' + str(1) + '\n')

    neg = os.listdir(negative_path)
    for i in range(len(neg)):
        if not re.match("(.*)npy", neg[i]): continue
        chrnum = int(neg[i].split('_')[0][3:])
        if chrnum % 2 == 0 or chrnum == 7 or chrnum == 19: continue
        dat = negative_path + neg[i]
        name = getname(dat)
        img = np.load(dat)
        if len(img) > 0:
            if type == 'indel':
                model = load_model('indel.m')
            else:
                model = load_model('snp.m')
            # y_pred = model.predict_classes(img, verbose=0).T[0]
            # sumn += len(y_pred)
            # for j in range(len(y_pred)):
            #     if y_pred[j] == 1:
            #         fp += 1
            #         f.write(name[j] + '\t' + str(0) + '\n')
            y_pred = model.predict_proba(img, verbose=0)
            sumn += len(y_pred)
            for j in range(len(y_pred)):
                if y_pred[j][1] + y_pred[j][2] > 0.5:
                    fp += 1
                    f.write(name[j] + '\t' + str(0) + '\n')
    f.close()
    print('Positive: ' + str(sump))
    print('Negative: ' + str(sumn))
    print('True positive rate: ' + str(float(tp) / sump))
    print('False positive: ' + str(fp))

def cnn():
    model = Sequential()
    model.add(Conv2D(32, 3, 3, activation='relu', border_mode='same', input_shape=(img_rows, img_cols, 3)))
    model.add(MaxPooling2D(pool_size=(2, 2), strides=(2, 2), border_mode='same'))
    model.add(Dropout(0.3))
    model.add(BatchNormalization())

    model.add(Conv2D(32, 3, 3, activation='relu', border_mode='same'))
    model.add(MaxPooling2D(pool_size=(2, 2), strides=(2, 2), border_mode='same'))
    model.add(Dropout(0.3))
    model.add(BatchNormalization())

    model.add(Conv2D(32, 3, 3, activation='relu', border_mode='same'))
    model.add(MaxPooling2D(pool_size=(2, 2), strides=(2, 2), border_mode='same'))
    model.add(Dropout(0.3))
    model.add(BatchNormalization())

    model.add(Conv2D(32, 3, 3, activation='relu', border_mode='same'))
    model.add(Dropout(0.3))
    model.add(BatchNormalization())

    model.add(Conv2D(32, 3, 3, activation='relu', border_mode='same'))
    model.add(Dropout(0.3))
    model.add(BatchNormalization())

    model.add(Conv2D(512, 1, 1, activation='relu', border_mode='same'))
    model.add(Flatten())
    # model.add(Dense(1, activation='sigmoid'))
    model.add(Dense(3, activation='softmax'))

    sgd = SGD(lr=0.005, decay=1e-6, momentum=0.9, nesterov=True)
    model.compile(loss='categorical_crossentropy', optimizer=sgd, metrics=['accuracy'])
    # model.compile(loss='binary_crossentropy', optimizer=sgd, metrics=['accuracy'])
    return model

def getlabel(y):
    yy = []
    for i in range(len(y)):
        if y[i][0] == 1:
            yy.append(0)
        elif y[i][1] == 1:
            yy.append(1)
        else:
            yy.append(2)
    return yy

def training(x_file, y_file, type):
    x_train = np.load(x_file)
    y_train = np.load(y_file)
    yy = getlabel(y_train)
    skf = StratifiedKFold(n_splits=5, random_state=1991, shuffle=True)
    for train_index, test_index in skf.split(x_train, yy):
        trn_x, val_x = x_train[train_index], x_train[test_index]
        trn_y, val_y = y_train[train_index], y_train[test_index]
        model = cnn()
        cal1 = ModelCheckpoint('best_weights.{epoch:02d}-{val_loss:.4f}.hdf5', monitor='val_loss',
                               save_best_only=True)
        cal2 = ReduceLROnPlateau(monitor='val_loss', factor=0.9, patience=2, min_lr=0.001)
        cal3 = EarlyStopping(monitor='val_loss', patience=5, mode='auto')
        model.fit(trn_x, trn_y, batch_size=200, nb_epoch=50, #callbacks=[cal1, cal2],
                  validation_split=0.0, shuffle=True, validation_data=(val_x, val_y))
        if type == 'indel':
            model.save('indel.m')
        else:
            model.save('snp.m')
        break


tye = 'indel'
print('Begin')
print(time.ctime())
training('x_train.npy', 'y_train.npy', tye)
print(time.ctime())
validate(tye)
print(time.ctime())
print('End')