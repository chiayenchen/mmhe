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
or download the entire repo from github website (`https://github.com/chiayenchen/mmhe`).


### How to use `mmhe` single GRM version
* Python version `mmhe.py`

  Usage: `python ./mmhe.py --grm my_grm_prefix --pheno my.pheno --mpheno 1 --covar my.covar`

  You can get the input files descrption with `./mmhe.py --help`.
  `mmhe.py` can take a pre-computed gentic relationship matrix (GRM), a phenotype file with multiple phenotypes, and a covaraite file (usually contains pricipal components for ancestry adjustment). Specifications of these files are listed below.

  1. The phenotype file follows the GCTA phenotype file format (same as PLINK phenotype format). The first 2 columns are the family ID and individual ID of the subjects included. These IDs are used to match the phenotype to the GRM. _Make sure these IDs correspond to the IDs in the genotype file used to calculate GRM._ IF you have multiple phenotypes in the file, you can specify which phenotype to use in the current analysis by `--mpheno`.

  2. The covariate file follows the GCTA `--qcovar` file format. The first 2 columns are the family ID and individual ID of the subjects included. These IDs are used to match the phenotype to the GRM. _Make sure these IDs correspond to the IDs in the genotype file used to calculate GRM._ Note that all covarites in the file will to be included in the analysis and all covariates are treated as continuous variables.

  3. The GRM follows GCTA binary GRM format (`PREFIX.grm.bin`, `PREFIX.grm.N.bin` and `PREFIX.grm.id`). However, `mmhe.py` only requires `PREFIX.grm.bin` and `PREFIX.grm.id`.

    Link to GCTA: `http://cnsgenomics.com/software/gcta/index.html`

    Link to PLINK2: `https://www.cog-genomics.org/plink2`

  The output of `mmhe.py` is the point esitamate and standard error of SNP-based heritabilty. The computation time is also provided.

* Matlab version `mmhe.m`

  Once read in the GRM, phenotype, and covariates in Matlab as matrices, `mmhe.m` can give the point esitamate and standard error of SNP-based heritabilty. Specifications of the input data format are listed below.

    1. y: n_subj x 1 vector of phenotype
    2. X: n_subj x n_cov matrix of covariates
    3. K: n_subj x n_subj matrix of empirical genetic similarity matrix
    (n_subj is the number of subjects anlyzed)

  The outputs from `mmhe.m` are

    1. h2: SNP heritability estimate
    2. se: standard error estimate of h2


### How to use `mmhe` vector stacking version
* Python version `mmhe_col.py`

  Usage: `python ./mmhe.py --grmdir my_grm_dir my_grm_prefix --pheno my.pheno --mpheno 1 --covar my.covar`

  For a dataset that has more than 100,000 subjects (or any *large* dataset), the `mmhe_col.py` can load the GRM by blockes.

  Specifications of the input data format are listed below.

    1. Phenotype file and covariate file follow the GCTA format (see description in the `mmhe` single GRM version).

    2. Block GRM files with `PREFIX.grm.id` file.

    The bock GRM files are n_subj x k matrices that are the columns of the *full* GRM matrix. Typically you would have `n_subj/k` block GRM files each with k columns of the full GRM and 1 block GRM files that will take less than `k` columns at the end of the GRM. These files should be plain text file and are saved as PREFIX.{block_num}.grm (e.g., PREFIX.1.grm, PREFIX.2.grm, ...) in 1 directory.

    The `PREFIX.grm.id` file shoud follow the format of GCTA GRM file format, with first column family ID and second column individual ID.

* Matlab version `mmhe_col.m`

  For a dataset that has more than 100,000 subjects (or any *large* dataset), the `mmhe_col.m` can load the GRM by blockes.

  Specifications of the input data format are listed below.
  
    1. y: n_subj x 1 vector of phenotype
    2. X: n_subj x n_cov matrix of covariates
    3. grm_dir: directory where block columns of the empirical genetic similarity matrix can be found; we have assumed here that each block column variable K is save as GRM_col{col_num}.mat (e.g., GRM_col1.mat, GRM_col2.mat, ...) in this directory.
    blk_size: size (number of columns) of each block
  
  The outputs from `mmhe_col.m` are
    
    1. h2: SNP heritability estimate
    2. se: standard error estimate of h2
    

### Support
Please contact Tian Ge (tge1@mgh.harvard.edu) or Chia-Yen Chen (chiayen.chen@mgh.harvard.edu) for any questions and comments.
