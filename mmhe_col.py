#!/usr/bin/env python
# =======================================
# matching-of-moment (MoM) method for heritability estimation
# input --
# y: n_subj x 1 vector of phenotype
# X: n_subj x n_cov matrix of covariates
# K: n_subj x n_subj matrix of empirical genetic similarity matrix
# output --
# h2: SNP heritability estimate
# se: standard error estimate of h2
# =======================================
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
from numpy import transpose as t
from numpy import dot as dot
from numpy import trace as tr
# import struct
import os
# import numpy as np
# import gzip

import time
start_time = time.time()

# =======================================
# parse command line input
INFO = 'Matching-of-moment (MoM) method for heritability estimation.\nE.g. mmhe.py\n--grmdir my_grm_dir my_grm_prefix block_size\n--pheno my.pheno\n--mpheno 1\n--covar my.covar\n\nThis version of mmhe will load block columns of the entire GRM from the specified directory.\n Each block column, which is a n_subj x block_size matrix, is save as PREFIX.{block_num}.grm (e.g., PREFIX.1.grm, PREFIX.2.grm, ...) in this directory.\n\nNote that you also need to provide a PREFIX.grm.id file to specify the subject ID in the block GRM files.\n\nAlso, all subjects must be included in the GRM, phenotype and covariate files.'
parser = argparse.ArgumentParser(description=INFO, formatter_class=RawTextHelpFormatter)

# parse PREFIX.grm.gz and PREFIX.grm.id file
# parser.add_argument('--grm-gz', type=str, help='Read GCTA format grm files, including PREFIX.grm.gz and PREFIX.grm.id [required]', required=True)

# parse PREFIX.grm.bin, PREFIX.grm.N.bin, and PREFIX.grm.id
parser.add_argument('--grmdir', nargs=2, help='Directory where block columns of the empirical genetic similarity matrix can be found.\n We have assumed here that each block column variable K is save as PREFIX.{block_num}.grm (e.g., PREFIX.1.grm, PREFIX.2.grm, ...) in this directory.\n The 2 arguments are\n1) the dirctory where blocks of GRM stored;\n2) PREFIX of the grm block files.\n\nNote that you also need to provide a PREFIX.grm.id file to specify the subject ID in the block GRM files.\n\nAlso, all subjects must be included in the GRM, phenotype and covariate files.', required=True)

# parse phenotype file
parser.add_argument('--pheno', type=str, help='Read PLINK format phenotype file [required]\nIf --mpheno is not specified then then 3rd column (the 1st phenotype) will be used.', required=True)
# specify the number of column for phenotype file
parser.add_argument('--mpheno', type=int, help='Specify which phenotype to use from phenotype file (1 phenotype only)')
parser.set_defaults(mpheno=1)
# parse covariate file
parser.add_argument('--covar', type=str, help='Read PLINK format covariate file.')
parser.set_defaults(covar="NULL")

args = parser.parse_args()
# print args

# =======================================
# provide a ID list for the GRM block file
# =======================================
# BinFileName = args.grm+".grm.bin"
IDFileName = args.grm[0]+"/"+args.grm[1]+".grm.id"
# NFileName = args.grm+".grm.N.bin"

id_list = []
with open(IDFileName, "r") as id:
    for i in id:
        itmp = i.rstrip().split()
        id_list.append(itmp[0]+":"+itmp[1])

n_subj = len(id_list)
nn = n_subj*(n_subj+1)/2

# =======================================
# read in coavriates
# =======================================
if args.covar == "NULL":
    X = np.ones(n_subj).reshape(n_subj, 1.0)
    n_cov = 1
else:
    covar_dic = {}
    # with open("./test.phen", "r") as X:
    with open(args.covar, "r") as X:
        for i in X:
            itmp = i.rstrip().split()
            covar_dic[itmp[0]+":"+itmp[1]] = [float(x) for x in itmp[2:]].insert(0, 1.0)
        n_cov = len(itmp) - 1

    covar_list = []
    for i in id_list:
        if i in covar_dic.keys():
            covar_list.append(covar_dic[i])

    X = np.array(covar_list).reshape(n_subj, n_cov)

# =======================================
# read in phenotype
# =======================================
pheno_dic = {}
# with open("./test.phen", "r") as y:
with open(args.pheno, "r") as y:
    for i in y:
        itmp = i.rstrip().split()
        pheno_dic[itmp[0]+":"+itmp[1]] = float(itmp[args.mpheno+1])

pheno_list = []
for i in id_list:
    if i in pheno_dic.keys():
        pheno_list.append(pheno_dic[i])

y = np.array(pheno_list).reshape(n_subj, 1)

# =======================================
# h2g
# =======================================
file_list = []
for file in os.listdir(args.grmdir[0]):
    if file.endswith('.grm') and file.startswith(args.grmdir[1]):
        file_list.append(file)

trK = 0.0
trKK = 0.0
yK = np.empty([1, n_subj], dtype=float)
XK = np.empty([n_cov, n_subj], dtype=float)
for i in file_list:
    Ktmp = []
    with open(args.grmdir[0]+"/"+i, "r") as Ktmp_file:
        for i in Ktmp_file:
            Ktmp.append(i.rstrip().split())
    block_size = len(Ktmp)/n_subj
    Ktmp = np.array(Ktmp).reshape(n_subj, block_size)
    for j in block_size:
        trK += Ktmp[j, j]

    trKK += np.sum(Ktmp*Ktmp)

    yK = np.append(yK, dot(t(y), Ktmp))
    XK = np.append(XK, dot(t(X), Ktmp), axis=1)

XX = dot(t(X), X)
XXinv = np.linalg.inv(XX)
Z = dot(X, XXinv)
yZ = dot(t(y), Z)
yPy = dot(t(y), y) - dot(dot(yZ, t(X)), y)

yZXK = dot(yZ, XK)
XKZ = dot(XK, Z)

yPKPy = dot(yK, y) - 2*dot(yZXK, y) + dot(dot(yZXK, X), t(yZ))
trPK = trK - tr(XKZ)
trPKPK = trKK - 2*tr(dot(dot(XXinv, XK), t(XK))) + tr(XKZ*XKZ)

S = np.array([trPKPK, trPK, trPK, n_subj-n_cov]).reshape(2, 2)
q = np.array([yPKPy, yPy]).reshape(2, 1)
Vc = dot(np.linalg.inv(S), q)

Vc[Vc < 0] = 0
s = trPKPK - trPK*trPK/(n_subj-n_cov)

h2 = np.asscalar(max(min(Vc[0]/sum(Vc), 1), 0))
se = pow(2/s, 0.5)
print h2, se

print("--- %s seconds ---" % (time.time() - start_time))
