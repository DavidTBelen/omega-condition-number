% This file contains some test matrices of sizes between 30000 
% and 6000 and discards the ones that are "trivial".
% We understand trivial if pcg with no preconditioner can solve the lsq
% problem in less than 10 seconds

%
clear all
seed=1;
rng(seed)

folder = fileparts(which(mfilename)); 
addpath(genpath(folder));

test_matrices = [ "sts4098","t2dal_e","bcsstk28", "msc04515","mhd4800b", "crystm01", "bcsstk16", "s3rmt3m3","ex15",'Muu',"bcsstk38",...
    'aft01', 'nd3k', 'fv1', 'fv3', 'bloweybq', 'bundle1', 'ted_B', 'msc10848', ...
    'bcsstk17', 't2dah_e', 'cbuckle', 'crystm02','Pres_Poisson', 'bcsstm25', 'Dubcova1',...
    'olafu','gyro', 'bodyy4', 'nd6k', 'bodyy5', 'bodyy6', 'raefsky4', 'LFAT5000'  , 'Trefethen_20000b',...
    't3dl_e', 'msc23052', 'crystm03', 'smt']; %, 'thread'];

small_matrices_list = [];



for ii=1:length(test_matrices)

    % Load data:
    load(test_matrices(ii))

    W = Problem.A;

    n = length(W); %size

    % Initialization:
    
    b = ones(n,1); %right hand side vector of ones
    tol = 1e-6;
    xinit = zeros(n,1);%starting point is the origin
    maxit = 100000;

    % No preconditioner:

    startW = tic;
    [xW,flagW,relresW,iterW,resvecW] = pcg(W,b,tol,maxit,[],[],xinit);
    stopW = toc(startW);
    disp('No preconditioner solver done!')

    if stopW > 10 || flagW >0
        small_matrices_list = [small_matrices_list, test_matrices(ii)];

    end




end





save('data/small_matrices_list.mat', 'small_matrices_list')