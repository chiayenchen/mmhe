function [h2, se] = mmhe_col(y, X, grm_dir, blk_size)
% memory efficient moment-matching method for SNP-based heritability estimation
%
% input --
% y: n_subj x 1 vector of phenotype
% X: n_subj x n_cov matrix of covariates
% grm_dir: directory where block columns of the empirical genetic
% similarity matrix can be found; we have assumed here that each block
% column variable K is save as GRM_col{col_num}.mat (e.g., GRM_col1.mat, GRM_col2.mat, ...) in this directory
% blk_size: size (number of columns) of each block
%
% output --
% h2: SNP heritability estimate
% se: standard error estimate of h2

n_subj = length(y);   % calculate the number of subjects n_subj
X = [ones(n_subj,1), X];   % add an intercept to the covariate matrix
n_cov = size(X,2);   % calculate the number of covariates n_cov

n_blk = ceil(n_subj/blk_size);   % calculate the total number of block columns

% iteratively load block columns of the empirical genetic similarity matrix
% to calculate trace(K), trace(K*K), y'*K, and X'*K 
trK = 0; trKK = 0; yK = zeros(1,n_subj); XK = zeros(n_cov,n_subj);
for ii = 1:n_blk
    load([grm_dir, '/GRM_blk', num2str(ii), '.mat']);   % load the ii-th block column
    if ii ~= n_blk
        trK = trK+trace(K((ii-1)*blk_size+1:ii*blk_size,:));
        trKK = trKK+sum(sum(K.^2));
        yK((ii-1)*blk_size+1:ii*blk_size) = y'*K;
        XK(:,(ii-1)*blk_size+1:ii*blk_size) = X'*K;
    else
        trK = trK+trace(K((ii-1)*blk_size+1:end,:));
        trKK = trKK+sum(sum(K.^2));
        yK((ii-1)*blk_size+1:end) = y'*K;
        XK(:,(ii-1)*blk_size+1:end) = X'*K;
    end
end

XX = X'*X;   % calculate X'*X
Z = X/(XX);   % calculate Z = X*(X'*X)^(-1)
yZ = y'*Z;   % calculate y'*Z
yPy = y'*y-yZ*X'*y;   % calculate y'*P*y, where P = I-X*(X'*X)^(-1)*X'

yZXK = yZ*XK;   % calculate y'*Z*X'*K
XKZ = XK*Z;   % calculate X'*K*Z

yPKPy = yK*y-2*yZXK*y+yZXK*X*yZ';   % calculate y'*P*K*P*y
trPK = trK-trace(XKZ);   % calculate trace(P*K)
trPKPK = trKK-2*trace(XX\XK*XK')+trace(XKZ*XKZ);   % calculate trace(P*K*P*K)

S = [trPKPK, trPK; trPK, n_subj-n_cov];   % S is a 2x2 matrix; first row: trace(P*K*P*K), trace(P*K); second row: trace(P*K), n_subj-n_cov
q = [yPKPy; yPy];   % q is a 2x1 vector
Vc = S\q;   % calculate variance component parameter estimates; Vc = S^(-1)*q

Vc(Vc<0) = 0;   % set neagtive variance component parameters to zero
s = trPKPK-trPK^2/(n_subj-n_cov);   % calculate s = trace(P*K*P*K)-[trace(P*K)]^2/(n_subj-n_cov)

h2 = max(min(Vc(1)/sum(Vc),1),0);   % calculate h2 estimate and restrict it to [0,1]
se = sqrt(2/s);   % standard error estimate