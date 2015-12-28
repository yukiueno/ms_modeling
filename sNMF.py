# coding: UTF-8
from numpy import *  # noqa
import scipy
from scipy.stats import uniform
import numpy
import random


def logsumexp(log_p):
    vec = log_p
    vec = scipy.misc.logsumexp(vec)
    return(log_p - scipy.misc.logsumexp(vec))


def update_g(g, X, p, c_0, D, Ct, rho):
    D = float(D)
    Ct = float(Ct)
    c_0 = float(c_0)
    xp = numpy.matrix([range(p.shape[0])] * X.shape[0], dtype=numpy.float)
    for i in range(X.shape[0]):
        for j in range(p.shape[0]):
            xp[i, j] = numpy.sum(X[i, :] * p[j, i, :])
    return((1 - rho) * g + rho * (c_0 / X.shape[0] + D / Ct * xp))


def update_h(h, a, b, c_0, D, Ct, rho):
    pre = ((1 - rho) * h)
    D = float(D)
    Ct = float(Ct)
    c_0 = float(c_0)
    pos = rho * (c_0 + D / Ct * numpy.sum(a / b, axis=1))
    ans = numpy.array(range(64), dtype=numpy.float)
    for i in range(64):
        ans[i] = pre[i] + pos[i]
    return(ans)


def update_a(X, p, a_0):
    sumation = numpy.array([range(X.shape[1])] * p.shape[0], dtype=numpy.float)
    for k in range(p.shape[2]):
        for i in range(p.shape[0]):
            sumation[i, k] = numpy.sum(X[:, k] * p[i, :, k])
    xp = sumation
    return(a_0 + xp)


def update_b(g, h, b_0):
    b = numpy.sum(g / h, axis=0)
    b = b_0 + b
    b = b.T
    return b


def calc_p(params):
    E_ln_theta = numpy.matrix([range(10)] * 64, dtype=numpy.float)
    for j in range(params['a'].shape[0]):
        for i in range(params['a'].shape[1]):
            E_ln_theta[j, i] = scipy.special.psi(params['a'][j, i]) - numpy.log(params['b'][j])
    E_ln_beta = numpy.matrix([range(64)] * 192, dtype=numpy.float)
    for i in range(params['g'].shape[1]):
        for j in range(params['g'].shape[0]):
            E_ln_beta[j, i] = scipy.special.psi(params['g'][j, i]) - numpy.log(params['h'][i])
    test = numpy.array([[range(E_ln_theta.shape[1])] * E_ln_beta.shape[0]] * 64, dtype=numpy.float)
    ppp = numpy.matrix([range(E_ln_theta.shape[1])] * E_ln_beta.shape[0], dtype=numpy.float)
    for i in range(E_ln_beta.shape[0]):
        for j in range(E_ln_theta.shape[1]):
            ooo = numpy.array(range(64), dtype=numpy.float)
            for k in range(64):
                ooo[k] = E_ln_beta[i, k] + E_ln_theta[k, j]
            ppp[i, j] = numpy.min(ooo)
    for i in range(E_ln_beta.shape[0]):
        for j in range(E_ln_theta.shape[1]):
            wa = numpy.array(range(64), dtype=numpy.float)
            for k in range(64):
                wa[k] = ((E_ln_beta[i, k] + E_ln_theta[k, j]) - ppp[i, j])
            test[:, i, j] = numpy.exp(wa - scipy.misc.logsumexp(wa))
    return(test)


def snmf(X, K, a_0=None, b_0=None, c_0=None,
         len_Ct=10, t_0=1, kappa=0.5, rho=1,
         n_itr=1e2, conv_l=0.1, conv_g=0.1):

    # init variational distribution
    numpy.random.seed()
    rv = uniform(loc=0.0, scale=1.0)
    g = numpy.matrix([rv.rvs(size=K)] * X.shape[0], dtype=numpy.float)
    for i in range(K):
        for j in range(X.shape[0]):
            g[j, i] = rv.rvs(size=1)[0]
    h = rv.rvs(size=K)
    a = numpy.matrix([rv.rvs(size=X.shape[1])] * K, dtype=numpy.float)
    for i in range(X.shape[1]):
        for j in range(K):
            a[j, i] = rv.rvs(size=1)[0]
    b = rv.rvs(size=K)
    numpy.c_[b]
    numpy.c_[h]

    # set parameters
    D_seq = X[1:, 0]
    params = dict(g=g, h=h, a=a, b=b)
    if a_0 is None:
        a_0 = 1.0 / K
    if b_0 is None:
        b_0 = 1.0 / K
    if c_0 is None:
        c_0 = 0.05 * X.shape[0]
    hp = dict(a_0=a_0, b_0=b_0, c_0=c_0)

    lb = numpy.array(range(int(n_itr)), dtype=numpy.float)

    for itr in (range(int(n_itr))):
        # selecting subset of the data
        Ct = random.sample(D_seq, len_Ct)
        X_sub = X[:, Ct]
        local_params = dict(a=params['a'][:, Ct], b=params['b'])

        for i in (range(int(1e3))):
            # compute p
            p = calc_p(dict(g=params['g'], h=params['h'], a=local_params['a'], b=local_params['b']))

            # def update_local variables
            local_params['a'] = update_a(X_sub, p, hp['a_0'])
            local_params['b'] = update_b(params['g'], params['h'], hp['b_0'])

            if i > 1:
                print "Θ更新中、現在の尤度差：%f" % numpy.sum(numpy.absolute(p - pre_p))
            if i > 1 and numpy.sum(numpy.absolute(p - pre_p)) < conv_l:
                break
            pre_p = p

        # def update_global variables
        params['g'] = update_g(params['g'], X_sub, p, hp['c_0'], X.shape[1], len_Ct, rho)
        params['h'] = update_h(params['h'], local_params['a'], local_params['b'], hp['c_0'], X.shape[1], len_Ct, rho)
        params['a'][:, Ct] = local_params['a']
        params['b'] = local_params['b']

        # def update_rho
        rho = (t_0 + itr) ** (-kappa)

        lb[itr] = numpy.sum(numpy.absolute(X - (params['g'] / params['h']) * (params['a'] / params['b'])))
        print "%d th iteration : %f" % (itr, lb[itr])
        if itr > 1 and numpy.absolute(lb[itr] - lb[itr - 1]) < conv_g:
            break

    return(dict(params=params, lb=lb))
