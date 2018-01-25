from keras.layers.normalization import BatchNormalization
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import MinMaxScaler
from sklearn.mixture import GaussianMixture
from sklearn.metrics import accuracy_score
from sklearn.decomposition import PCA
from sklearn.externals import joblib
from sklearn.neighbors import KDTree
from keras.layers.core import Dense
from keras.models import Sequential
from keras.models import load_model
from sklearn.cluster import KMeans
from keras.optimizers import SGD
from sklearn import svm
# import xgboost as xgb
import numpy as np
import random

datapath = '/home/linfengye/mutations/data/data_new/'
snp_p = datapath + 'bed.na878.mut.snp.nofilt.TP'
snp_n = datapath + 'bed.na878.mut.snp.nofilt.FP'
indel_p = datapath + 'bed.na878.mut.indel.nofilt.TP'
indel_n = datapath + 'bed.na878.mut.indel.nofilt.FP'

def upsample(x, n):
    x_new = []
    tree = KDTree(x, leaf_size=2)
    for i in range(len(x)):
        x_new.append(x[i])
        tt = np.array(x[i])
        dist, index = tree.query([tt], k=n)
        for j in range(index.shape[0]):
            ttt = np.array(x[index[0][j]])
            tttt = tt + (ttt - tt) * np.random.random(size=14)
            x_new.append(list(tttt))
    return x_new

def chara():
    x = list(np.load('x_train.npy'))
    x1 = list(np.load('x_validate.npy'))
    x.extend(x1)
    y = list(np.load('y_train.npy'))
    y1 = list(np.load('y_validate.npy'))
    y.extend(y1)
    y = np.array(y)
    pca = PCA(n_components=3)
    newx = pca.fit_transform(x)
    newx = np.array(newx)
    min_max_scaler = MinMaxScaler()
    newx = min_max_scaler.fit_transform(newx)
    newx = np.array(newx)
    np.save('x.npy', newx)
    np.save('y.npy', y)

def net():
    model = Sequential()
    model.add(BatchNormalization(input_shape=(14,)))
    model.add(Dense(8, activation='relu'))
    model.add(Dense(8, activation='relu'))
    model.add(Dense(4, activation='relu'))
    model.add(Dense(1, activation='sigmoid'))
    sgd = SGD(lr=0.1, decay=1e-6, momentum=0.9, nesterov=True)
    model.compile(loss='binary_crossentropy', optimizer=sgd, metrics=['accuracy'])
    return model

def visualization():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    x = np.load('x.npy')
    y = np.load('y.npy')
    ax = plt.subplot(111, projection='3d')
    for i in range(len(x)):
        if y[i] == 1:
            color = 'r'
        else:
            color = 'g'
        ax.scatter(x[i][0], x[i][1], x[i][2], c=color)
    ax.set_zlabel('Z')
    ax.set_ylabel('Y')
    ax.set_xlabel('X')
    plt.savefig("example.tiff")

def data(file, type):
    inf = []
    arr = []
    sf = open(file, 'r')
    while True:
        line = sf.readline()
        if not line: break
        stri = line.split()
        chr = stri[0]
        if chr == 'Chr' or float(stri[-1]) < -3 or float(stri[-2]) > -8: continue
        if type == 'indel':
            if len(stri[3]) == 1: continue
        else:
            if len(stri[3]) > 1: continue
        poi = stri[1]
        out = chr + '_' + poi
        inf.append(out)
        cstr = stri[4:]
        for i in range(len(cstr)):
            if cstr[i] == 'N':
                cstr[i] = 0
            cstr[i] = float(cstr[i])
        arr.append(cstr)
    arr = np.array(arr)
    return inf, arr

def validate(type, id, file):
    clf = joblib.load(type+'.m')
    inf, arr = data(file, type)
    y_pred = clf.predict(arr)
    f = open(id+'.'+type, 'w')
    for i in range(len(y_pred)):
        if y_pred[i] == 1:
            f.write(inf[i] + '\n')
    f.close()

###########################################################################

def getSample(file):
    x_train = []
    x_validate = []
    inf = []
    f = open(file)
    while True:
        line = f.readline()
        if not line: break
        cstr = line.split()
        chr = cstr[0][3:]
        if chr == '19': continue
        cp = chr + '_' + cstr[1]
        chr = int(chr)
        # if float(cstr[-1]) < -3 or float(cstr[-2]) > -8: continue	# filt for snp
        cstr = cstr[4:]
        for i in range(len(cstr)):
            if cstr[i] == 'N':
                cstr[i] = 0
            cstr[i] = float(cstr[i])
        if chr < 19:
            x_train.append(cstr)
        else:
            x_validate.append(cstr)
            inf.append(cp)
    return x_train, x_validate, inf

def getTrain(train_p, train_n):
    num_p = len(train_p)
    num_n = len(train_n)
    poi_p = 0
    poi_n = 0
    label = 1
    train_x = []
    train_y = []
    while label:
        if random.random() < 0.5:
            train_x.append(train_n[poi_n])
            train_y.append(0)
            poi_n += 1
        else:
            train_x.append(train_p[poi_p])
            train_y.append(1)
            poi_p += 1
        if poi_n == num_n:
            for i in range(poi_p, num_p):
                train_x.append(train_p[i])
            train_y.extend([1] * (num_p - poi_p))
            label = 0
        if poi_p == num_p:
            for i in range(poi_n, num_n):
                train_x.append(train_n[i])
            train_y.extend([0] * (num_n - poi_n))
            label = 0
    train_x = np.array(train_x)
    train_y = np.array(train_y)
    return train_x, train_y

def sample(type):
    if type  == 'snp':
        train_p, validate_p, inf_p = getSample(snp_p)
        train_n, validate_n, inf_n = getSample(snp_n)
    else:
        train_p, validate_p, inf_p = getSample(indel_p)
        train_n, validate_n, inf_n = getSample(indel_n)
    inf = inf_p + inf_n
    x_train, y_train = getTrain(train_p, train_n)
    x_validate = validate_p + validate_n
    y_validate = len(validate_p) * [1] + len(validate_n) * [0]
    x_validate = np.array(x_validate)
    y_validate = np.array(y_validate)
    print(x_train.shape)
    print(x_validate.shape)
    np.save('x_train.npy', x_train)
    np.save('y_train.npy', y_train)
    np.save('x_validate.npy', x_validate)
    np.save('y_validate.npy', y_validate)
    f = open('validate.list', 'w')
    for i in range(len(inf)):
        f.write(inf[i] + '\n')

def train(type, is_train):
    x_train = np.load('x_train.npy')
    y_train = np.load('y_train.npy')
    x_validate = np.load('x_validate.npy')
    y_validate = np.load('y_validate.npy')
    # model = net()
    # model.fit(x_train, y_train, batch_size=200, nb_epoch=6, validation_split=0.0, shuffle=True)
    # y_pred = model.predict_classes(x_validate, verbose=0)
    # clf = svm.SVC(C=1.0, kernel='rbf', degree=3, gamma='auto', coef0=0.0, cache_size=51200)
    # clf = GradientBoostingClassifier(learning_rate=0.1, n_estimators=100)
    # clf = xgb.XGBRegressor(max_depth=6, n_estimators=1500, min_child_weight=65, learning_rate=0.036,
    #                        subsample=0.90, colsample_bytree=0.90, seed=1991)
    # clf = KNeighborsClassifier(n_neighbors=100, algorithm='auto', leaf_size=30, n_jobs=20)
    # clf = GaussianMixture(n_components=2, max_iter=200, n_init=5, random_state=1991)
    # clf = KMeans(n_clusters=2, max_iter=200, n_init=1, random_state=1991)
    # clf.fit(x_train)
    cutoff = 0.5
    if is_train:
        clf = RandomForestClassifier(n_estimators=600, max_features=None, min_samples_split=20,
                                     min_samples_leaf=50, oob_score=True, n_jobs=20, class_weight=None)
        clf.fit(x_train, y_train)
        joblib.dump(clf, type+'.m')
    else:
        clf = joblib.load(type + '.m')
    y_pred = clf.predict_proba(x_validate)
    f = open('validate.list', 'r')
    inf = []
    while True:
        line = f.readline()
        if not line: break
        inf.append(line.split()[0])
    fr = open(type + '_result.xls', 'w')
    fw = open(type + '_false.xls', 'w')
    p = 0; n = 0
    tp = 0; fp = 0
    for i in range(len(y_pred)):
        fr.write(str(y_pred[i][1]) + '\t' + str(y_validate[i]) + '\n')
        if y_validate[i] == 1:
            p += 1
        else:
            n += 1
        if (y_pred[i][1] - cutoff) * (y_validate[i] - cutoff) >= 0 and y_validate[i] == 1:
            tp += 1
        elif (y_pred[i][1] - cutoff) * (y_validate[i] - cutoff) < 0 and y_validate[i] == 0:
            fp += 1
            fw.write(inf[i] + '\t' + str(y_validate[i]) + '\n')
        elif (y_pred[i][1] - cutoff) * (y_validate[i] - cutoff) < 0 and y_validate[i] == 1:
            fw.write(inf[i] + '\t' + str(y_validate[i]) + '\n')
    fr.close()
    fw.close()
    print('Training: ' + str(len(x_train)) + ', validate: ' + str(len(x_validate)))
    print('Positive: ' + str(p) + ', negative: ' + str(n))
    print('True positive: ' + str(tp))
    print('False positive: ' + str(fp))
    print('True positive rate: ' + str(float(tp)/p))


type = 'indel'
is_train = 1
if is_train:
    sample(type)
train(type, is_train)