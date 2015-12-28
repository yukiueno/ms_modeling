# coding: UTF-8
from numpy import *


def difcost(a, b):
    dif = 0
    # 行列の全ての行と列をループする。
    for i in range(shape(a)[0]):
        for j in range(shape(a)[1]):
            # 差を足し合わせる
            dif += pow(a[i, j] - b[i, j], 2)
    return dif


def factorize(v, pc=10, iter=50):
    ic = shape(v)[0]  # shape関数は行列の行と列の数を返す。
    fc = shape(v)[1]

    # 重みと特徴の行列をランダムな値で初期化
    w = matrix([[random.random() for j in range(pc)] for i in range(ic)])
    h = matrix([[random.random() for i in range(fc)] for i in range(pc)])

    # 最大でiterの回数だけ操作を繰り返す:
    for i in range(iter):
        wh = w * h

        # 現在の差を計算
        cost = difcost(v, wh)

        # 行列が完全に因子分解されたら終了
        if cost == 0:
            break

        # 特徴の行列を更新
        hn = (transpose(w) * v)
        hd = (transpose(w) * w * h)

        h = matrix(array(h) * array(hn) / array(hd))

        # 重みの行列を更新
        wn = (v * transpose(h))
        wd = (w * h * transpose(h))

        w = matrix(array(w) * array(wn) / array(wd))

    return w, h
