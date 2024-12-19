% This file compares different preconditioners for solving linear systems
% generated with some set of test matrices
% We use the pcg as a solver

clear 

seed=1;
rng(seed);

folder = fileparts(which(mfilename)); 
addpath(genpath(folder));

% Choose the data set: 
which_experiment = 'small_matrices';
data = load([which_experiment,'_list']); %array with name of matrices
problems = data.small_matrices_list;



nproblems = length(problems);

% Arrays for storing the output:

output.nproblems = nproblems;

output.problem_sizes = zeros(nproblems,1);
output.problem_nnz = zeros(nproblems,1);
output.problem_density = zeros(nproblems,1);


% No preconditioner:
output.no_prec.it = zeros(nproblems,1);
output.no_prec.time = zeros(nproblems,1);
output.no_prec.relres = zeros(nproblems,1);
output.no_prec.res = zeros(nproblems,1);
output.no_prec.fails = zeros(nproblems,1);


% Diagonal preconditioner:
output.diag_prec.it = zeros(nproblems,1);
output.diag_prec.time = zeros(nproblems,1);
output.diag_prec.time_prec = zeros(nproblems,1);
output.diag_prec.relres = zeros(nproblems,1);
output.diag_prec.res = zeros(nproblems,1);
output.diag_prec.fails = zeros(nproblems,1);


% Diagonal + block triangular preconditioner:
output.dbt_prec.it = zeros(nproblems,1);
output.dbt_prec.time = zeros(nproblems,1);
output.dbt_prec.time_prec = zeros(nproblems,1);
output.dbt_prec.relres = zeros(nproblems,1);
output.dbt_prec.res = zeros(nproblems,1);
output.dbt_prec.fails = zeros(nproblems,1);

% Incomplete Cholesky preconditioner:
output.ichol_prec.it = zeros(nproblems,1);
output.ichol_prec.time = zeros(nproblems,1);
output.ichol_prec.time_prec = zeros(nproblems,1);
output.ichol_prec.relres = zeros(nproblems,1);
output.ichol_prec.res = zeros(nproblems,1);
output.ichol_prec.fails = zeros(nproblems,1);

% Incomplete Cholesky preconditioner 2:
output.ichol_2_prec.it = zeros(nproblems,1);
output.ichol_2_prec.time = zeros(nproblems,1);
output.ichol_2_prec.time_prec = zeros(nproblems,1);
output.ichol_2_prec.relres = zeros(nproblems,1);
output.ichol_2_prec.res = zeros(nproblems,1);
output.ichol_2_prec.fails = zeros(nproblems,1);


for pp = 1:nproblems
    
    fprintf('#########  Problem %i  ######### \n', pp)
    % Load data:
    load(problems(pp))

    W = Problem.A;


    n = length(W); %size

    nnzW = nnz(W);

    output.problem_sizes(pp) = n;
    output.problem_nnz(pp) = nnzW;
    output.problem_density(pp) = nnzW/n^2;


    k =  ceil((.5+sqrt(1+4*2*nnzW*0.1)*.5))+1; %dimension for the triangular block of the preconditioner
    
    % Compute optimal preconditioners:


    % Diagonal preconditioner:
    start_diag_prec = tic;
    D = diag_prec(W);
    stop_dp = toc(start_diag_prec);
    %disp(' Diag_prec done! \n')
    WD = D'*W*D;
    WD = (WD+WD')/2;
    output.diag_prec.time_prec(pp) = stop_dp;

    

    % Upper-block-triangular preconditioner:
    start_bt_prec = tic;
    U = i_upper_tri_preconditioner(W,k);
    stop_btp = toc(start_bt_prec);
    %disp('Solver: block_trir_prec done! \n')
    WU = U'*W*U;
    WU = (WU+WU')/2;
    output.dbt_prec.time_prec(pp) = stop_btp;



    % Incomplete Cholesky preconditioner
    start_ichol_prec = tic;
    alpha = max(sum(abs(W),2)./diag(W))-2;
    L = ichol(W, struct('type','ict','droptol',1e-3,'diagcomp',alpha));
    %H = inv(L');    
    stop_ichol_prec = toc(start_ichol_prec);
    output.ichol_prec.time_prec(pp) = stop_ichol_prec;

    % Incomplete Cholesky preconditioner (2)
    start_ichol_2_prec = tic;
    alpha = max(sum(abs(W),2)./diag(W))-2;
    alpha = alpha*0.01;
    L2 = ichol(W, struct('type','ict','droptol',1e-3,'diagcomp',alpha));
    %H = inv(L');    
    stop_ichol_2_prec = toc(start_ichol_2_prec);
    output.ichol_2_prec.time_prec(pp) = stop_ichol_2_prec;



    % Run solvers:

    % Initialization:
    
    b = ones(n,1); %right hand side vector of ones
    tol = 1e-6;
    xinit = zeros(n,1);%starting point is the origin
    maxit = 100000;

    % No preconditioner:

    startW = tic;
    [xW,flagW,relresW,iterW,resvecW] = pcg(W,b,tol,maxit,[],[],xinit);
    stopW = toc(startW);
    disp('1. No preconditioner solver done!')

    output.no_prec.it(pp) = iterW;
    output.no_prec.time(pp) = stopW;
    output.no_prec.relres(pp) = relresW;
    output.no_prec.res(pp) = norm(W*xW-b);
    if flagW > 0 || relresW > tol
        output.no_prec.fails(pp) = 1;
    end

    
    % Diagonal preconditioner:

    startD = tic;
    [yD,flagD,relresD,iterD,resvecD] = pcg(WD,D'*b,tol,maxit,[],[],xinit);
    xD = D*yD;
    stopD = toc(startD);
    disp('2. Diagonal preconditioner solver done!')

    output.diag_prec.it(pp) = iterD;
    output.diag_prec.time(pp) = stopD + stop_dp;
    output.diag_prec.relres(pp) = relresD;
    output.diag_prec.res(pp) = norm(W*xD-b);
    if flagD > 0 || relresD > tol
        output.diag_prec.fails(pp) = 1;
    end

    % Diagonal + block triangular preconditioner

    startU = tic;
    [yU,flagU,relresU,iterU,resvecU] = pcg(WU,U'*b,tol,maxit,[],[],xinit);
    xU = U*yU;
    stopU = toc(startU);
    disp('3. Diagonal + block triangular solver done!')


    output.dbt_prec.it(pp) = iterU;
    output.dbt_prec.time(pp) = stopU + stop_btp;
    output.dbt_prec.relres(pp) = relresU;
    output.dbt_prec.res(pp) = norm(W*xU-b);
    if flagU > 0 || relresU > tol
        output.dbt_prec.fails(pp) = 1;
    end

    % Incomplete Cholesky
    startC = tic;
    [xC,flagC,relresC,iterC,resvecC] = pcg(W,b,tol,maxit,L,L',xinit);
    stopC = toc(startC);
    disp('4. Incomplete cholesky solver done!')
    
    output.ichol_prec.it(pp) = iterC;
    output.ichol_prec.time(pp) = stopC + stop_ichol_prec;
    output.ichol_prec.relres(pp) = relresC;
    output.ichol_prec.res(pp) = norm(W*xC-b);
    if flagC > 0 || relresC > tol
        output.ichol_prec.fails(pp) = 1;
    end

    
    % Incomplete Cholesky 2
    startC_2 = tic;
    [xC_2,flagC_2,relresC_2,iterC_2,resvecC_2] = pcg(W,b,tol,maxit,L2,L2',xinit);
    stopC_2 = toc(startC_2);
    disp('5. Incomplete cholesky solver done!')
    
    output.ichol_2_prec.it(pp) = iterC_2;
    output.ichol_2_prec.time(pp) = stopC_2 + stop_ichol_2_prec;
    output.ichol_2_prec.relres(pp) = relresC_2;
    output.ichol_2_prec.res(pp) = norm(W*xC_2-b);
    if flagC_2 > 0 || relresC_2 > tol
        output.ichol_2_prec.fails(pp) = 1;
    end

    fprintf('\n')

    

end % end for pp

 filename =  "data/output_exp_" + string(which_experiment);
 save(filename,"output")
     

