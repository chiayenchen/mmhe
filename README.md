# mmhe
## Moment matching method for SNP-based heritability estimation

`mmhe` is an implemtation of moment-matching method for SNP-based heritability eitmation.


### Getting started
For the python scripts, you will need to install python and the packages required, including `argparse`, `numpy`, `os`, and `struct`.

For the Matlab scripts, you will need to install Matlab.

You can get `mmhe`by simply clone this repo with 
```  
git clone https://github.com/chiayenchen/mmhe.git
```
or download the entire repo from github website.


### How to use `mmhe` single GRM version
* Python version `mmhe.py`
You can get the input files descrption with `./mmhe.py --help`.
`mmhe.py` can take a pre-computed gentic relationship matrix (GRM), a phenotype file with multiple phenotypes, and a covaraite file (usually contains pricipal components for ancestry adjustment). Specifications of these files are listed below.

1. The GRM follows GCTA binary GRM format (`PREFIX.grm.bin`, `PREFIX.grm.N.bin` and `PREFIX.grm.id`). However, `mmhe.py` only requires `PREFIX.grm.bin` and `PREFIX.grm.id`. We recommend using `PLINK2` for calculating GRM. 
Link to GCTA: `http://cnsgenomics.com/software/gcta/index.html`
Link to PLINK2: `https://www.cog-genomics.org/plink2`

2. The phenotype file follows the GCTA phenotype file format (same as PLINK phenotype format). The first 2 columns are the family ID and individual ID of the subjects included. These IDs are used to match the phenotype to the GRM. _Make sure these IDs correspond to the IDs in the genotype file used to calculate GRM._ IF you have multiple phenotypes in the file, you can specify which phenotype to use in the current analysis by `--mpheno`.

3. The covariate file follows the GCTA `--qcovar` file format. The first 2 columns are the family ID and individual ID of the subjects included. These IDs are used to match the phenotype to the GRM. _Make sure these IDs correspond to the IDs in the genotype file used to calculate GRM._ Note that all covarites in the file will to be included in the analysis and all covariates are treated as continuous variables.

The output of `mmhe.py` is the point esitamate and standard error of SNP-based heritabilty. The computation time is also provided.

* Matlab version `mmhe.m`
Once read in the GRM, phenotype, and covariates in Matlab as matrices, `mmhe.m` can give the point esitamate and standard error of SNP-based heritabilty. Specifications of the input data format are listed below.

y: n_subj x 1 vector of phenotype
X: n_subj x n_cov matrix of covariates
K: n_subj x n_subj matrix of empirical genetic similarity matrix
(n_subj is the number of subjects anlyzed)

The outputs from `mmhe.m` are

h2: SNP heritability estimate
se: standard error estimate of h2

### How to use `mmhe` vector stacking version
* Python version `mmhe_col.py`


* Matlab version `mmhe_col.m`
For a dataset that has more than 100,000 subjects (or any *large* dataset), the `mmhe_col.m` can load the GRM by blockes to make the I/O more efficient.
Specifications of the input data format are listed below.

y: n_subj x 1 vector of phenotype
X: n_subj x n_cov matrix of covariates
grm_dir: directory where block columns of the empirical genetic
similarity matrix can be found; we have assumed here that each block
column variable K is save as GRM_col{col_num}.mat (e.g., GRM_col1.mat, GRM_col2.mat, ...) in this directory.
blk_size: size (number of columns) of each block

The outputs from `mmhe_col.m` are

h2: SNP heritability estimate
se: standard error estimate of h2


### Support
Please contact Tian Ge (tge1@mgh.harvard.edu) or Chia-Yen Chen (chiayen.chen@mgh.harvard.edu) for any questions and comments.