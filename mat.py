# coding: UTF-8
import sys
import read
# import NMF
import sNMF
from numpy import *  # noqa

M = read.read_mat_file(sys.argv)
# w, h = NMF.factorize(M, pc=4, iter=100)
w = sNMF.snmf(M, 64, a_0=None, b_0=None, c_0=None, len_Ct=10, t_0=1, kappa=0.5, rho=1, n_itr=1e2, conv_l=1e-2, conv_g=1e-2)

print w

"""
# 分解された行列を出力
print 'w='
for i in range(shape(w)[0]):
    for j in range(shape(w)[1]):
        print ('%3f' % w[i, j]),
    print
print
print 'h='
for i in range(shape(h)[0]):
    for j in range(shape(h)[1]):
        print ('%3f' % h[i, j]),
    print
print
"""
# 元々の行列を出力
print 'answer='
print M
print

# 分解された行列２つを掛け合わせた答えを出力
print 'w*h='
wh = w * h
for i in range(shape(wh)[0]):
    for j in range(shape(wh)[1]):
        print ('%3f' % wh[i, j]),
    print
