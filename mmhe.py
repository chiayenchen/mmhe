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
import struct
# import numpy as np
# import gzip

# =======================================
# parse command line input
INFO = 'Matching-of-moment (MoM) method for heritability estimation.\nE.g. mmhe.py\n--grm my_grm_prefix\n--pheno my.pheno\n--mpheno 1\n--covar my.covar\n\nAll files are assumed to be white space delimited.'
parser = argparse.ArgumentParser(description=INFO, formatter_class=RawTextHelpFormatter)

# parse PREFIX.grm.gz and PREFIX.grm.id file
# parser.add_argument('--grm-gz', type=str, help='Read GCTA format grm files, including PREFIX.grm.gz and PREFIX.grm.id [required]', required=True)

# parse PREFIX.grm.bin, PREFIX.grm.N.bin, and PREFIX.grm.id
parser.add_argument('--grm', type=str, help='Read GCTA binary format grm files, including PREFIX.grm.bin, PREFIX.grm.N.bin, and PREFIX.grm.id [required]', required=True)

# parse phenotype file
parser.add_argument('--pheno', type=str, help='Read PLINK format phenotype file. If --pheno-number is not specified then then 3rd column will be used [required]', required=True)
# specify the number of column for phenotype file
parser.add_argument('--mpheno', type=int, help='Specify which column to use for phenotype file (1 phenotype only)')
parser.set_defaults(mpheno=1)
# parse covariate file
parser.add_argument('--covar', type=str, help='Read PLINK format covariate file.', required=True)
parser.set_defaults(covar="NULL")

# estimate h2g
# parser.add_argument('--h2g', action='store_true', help='Estimate heritability')
# parser.set_defaults(h2g=False)

# specify output result file name
# parser.add_argument('--out', type=str, help='Specify output file name', required=True)

args = parser.parse_args()
# print args

# =======================================
# read in grm
# =======================================
BinFileName = args.grm+".grm.bin"
IDFileName = args.grm+".grm.id"
# NFileName = args.grm+".grm.N.bin"

# fileprefix="./test"
# BinFileName=fileprefix+".grm.bin"
# NFileName=fileprefix+".grm.N.bin"
# IDFileName=fileprefix+".grm.id"

id_list = []
with open(IDFileName, "r") as id:
    for i in id:
        itmp = i.rstrip().split()
        id_list.append(itmp[0]+":"+itmp[1])

n_subj = len(id_list)
nn = n_subj*(n_subj+1)/2

with open(BinFileName, "rb") as BinFile:
    BinFileContent = BinFile.read()

K = np.zeros((n_subj, n_subj))
grm_vals = list(struct.unpack("f"*nn, BinFileContent))
# grm_vals = [ind1*ind1, ind2*ind1, ind2*ind2, ind3*ind1, ind3*ind2, ind3*ind3, ...]
inds = np.tril_indices_from(K)
K[inds] = grm_vals
K[(inds[1], inds[0])] = grm_vals

# =======================================
# read in coavriates
# =======================================
if args.covar == "NULL":
    X = np.ones(n_subj)
else:
    covar_dic = {}
    # with open("./test.phen", "r") as X:
    with open(args.covar, "r") as X:
        for i in X:
            itmp = i.rstrip().split()
            covar_dic[itmp[0]+":"+itmp[1]] = [float(x) for x in itmp[2:]]
        n_cov = len(itmp) - 2

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
trK = np.trace(K)
trKK = np.sum(K*K)
yK = np.dot(np.transpose(y), K)
XK = np.dot(np.transpose(X), K)

XX = np.dot(np.transpose(X), X)
Z = X/XX
yZ = np.dot(np.transpose(y), Z)
yPy = np.dot(np.transpose(y), y) - np.dot(np.dot(yZ, np.transpose(X)), y)

yZXK = np.dot(yZ, XK)
XKZ = np.dot(XK, Z)

yPKPy = np.dot(yK, y) - 2*np.dot(yZXK, y) + np.dot(np.dot(yZXK, X), np.transpose(yZ))
trPK = trK - np.trace(XKZ)
trPKPK = trKK - 2*np.trace(np.dot(np.dot(np.linalg.inv(XX), XK), np.transpose(XK))) + np.trace(XKZ*XKZ)

S = np.array([trPKPK, trPK, trPK, n_subj-n_cov]).reshape(2, 2)
q = np.array([yPKPy, yPy]).reshape(2, 1)
Vc = np.dot(np.linalg.inv(S), q)

Vc[Vc < 0] = 0
s = trPKPK - trPK*trPK/(n_subj-n_cov)

h2 = np.asscalar(max(min(Vc[0]/sum(Vc), 1), 0))
se = pow(2/s, 0.5)
print h2, se
