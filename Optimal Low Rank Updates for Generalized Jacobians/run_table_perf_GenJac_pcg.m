%%% This file generates the table and performance profile in Section 4.2.
%%% 
clear
seed = 2;
rng(seed);
profile clear
profile off

options.savetable = true; %true for save table
options.saveplot = true; %true for sace plot


% sizesn = [2000];  %dimensions for tests
sizesn = [1000,2000,3000,5000];  %dimensions for tests
resultsn = zeros(length(sizesn),5,6);
optionsplot.sizesn = sizesn;

solver = 'pcg'; % 'lsqr' or 'cgs'
optionsplot.solver =  solver;

timegammastarn = zeros(length(sizesn),4);
totaltimegammastarn = zeros(length(sizesn),3);

%%% table name and performance profile figure name
tablename = 'tablejacob';   % name used for latex file as of Oct24/2024
perfprofnameiter = 'perfprofiter';  % name used for latex file as of aug15/23
perfprofnametime = 'perfprofomega';  % name used for latex file as of aug15/23
%  ['../../../condnumbandQPlatexfiles.d/perfprofiter', solver, '.png'];


np = 10; %number of problems
%np = 5; %number of problems
% np = 2; %number of problems

setting = zeros(length(sizesn),np,4); %for storing n,r and t
results = zeros(length(sizesn),np,5,6);  %for storing iterations, cond, residuals and time
timegammastar = zeros(length(sizesn),np,4); %for storing the time for computing gammastar
counterconst = zeros(length(sizesn),np); %for counting the gammastar outside [0,1]


for nn =1:length(sizesn)  %for dimension
    n = sizesn(nn);
    
    tol = 1e-12;
                    
    densityU = 1/log(n)*rand;
    densityA = .5/log(n)*rand;  % smaller than  for U as it gets squared
    for npp=1:np      % number of problems each dimension
        epsilon = 1e-7*(rand+1e-2);  % for regularizing A
        r = min(n-1,round(1+n*max(.5,rand))); % rank of A
        t = max(round(min(r/2,round(1+n*rand/2))),2);   % rank of update
        
        % storing values of n, r and t
        setting(nn,npp,1)=n;
        setting(nn,npp,2)=r;
        setting(nn,npp,3)=t;
        %
        A = sprandn(n,r,densityA); 
        U = sprandn(n,t,densityU);
        U(:,norms(U) == 0)= [];
        t = size(U,2);
        nus = full((norms(U).^2)');

        A = A*A';   
        A = (A+A')/2; % final A, symmetric and psd (singular)
        B = U*randn(t,1) +A(:,1:n)*randn(n,1);   % solve with same RHS for three systems
        %keyboard;

        startgamma=tic;   % for spectral method
        [q,Ddiag] = eig(full(A));
        Ds = sqrt(Ddiag+epsilon*speye(n));
        W = (q/Ds)'*U;
        
        nws = full((norms(W).^2)');

        invnus = nus.^(-1);
        Kinv = 1/n*diag(invnus) + repmat(invnus,1,t)/(n*(n-t)); % Sherwin-Morr.
        b0 = (trace(A)+n*epsilon) - (n*nus./nws);
        gammastar = Kinv*b0; %optimal gamma
        %%% Projecting to obtain the [0,1] 
        if max(gammastar)>1 || min(gammastar)<0
            fprintf('Gamma star outside [0,1] \n')
            counterconst(nn,npp) = 1;
            gammastar(gammastar>1) = 1;
            gammastar(gammastar<0) = 0;
        end
        stopgamma = toc(startgamma);
        timegammastar(nn,npp,1) = stopgamma;
        
        %Computing gammastar with Cholesky
        startgammaC = tic;
        tic
        R = chol(A+epsilon*speye(n));
        Wc = (R')\U;
        nwsc = full((norms(Wc).^2)');
        invnus = nus.^(-1);
        Kcinv = 1/n*diag(invnus) + repmat(invnus,1,t)/(n*(n-t)); % Sherwin-Morr.
        bc0 = (trace(A)+n*epsilon) - (n*nus./nwsc);
        gammastarC = Kcinv*bc0; %optimal gamma
        stopgammaC = toc(startgammaC);
        timegammastar(nn,npp,2) = stopgammaC;

        %Computing gammastar with approximation
        startgammaA = tic;
        nun = full((norms(U).^2)');
        gammastarA = trace(A)./((n-t)*nun);
        gammastarA(gammastarA>1) = 1;
        stopgammaA = toc(startgammaA);
        timegammastar(nn,npp,3) = stopgammaA;

        Aeps = A + (max(0,-min(eig(A)))+epsilon)*speye(n);
        A0 = Aeps;
        A1 = Aeps + U*U';

        %Computing gammastar = 1/u.^2
        startcompn = tic;
        nun = full((norms(U).^2)');
        gamman = 1./nun;
        gamman(gamman>1) = 1;
        stopcompn = toc(startcompn);
        timegammastar(nn,npp,4) = stopcompn;
        An = Aeps + U*diag(gamman)*U';
        An = (An+An')/2;
        As = Aeps + U*diag(gammastar)*U';
        As = (As+As')/2;
        Aa = Aeps + U*diag(gammastarA)*U';
        Aa = (Aa+Aa')/2;
        
        % storing the condest condition number
        results(nn,npp,1,1) = condest(A0);
        results(nn,npp,2,1) = condest(A1);
        results(nn,npp,3,1) = condest(An);
        results(nn,npp,4,1) = condest(As);
        results(nn,npp,5,1) = condest(Aa);
        % storing the omega condition number
        eig0 = eig(A0);
        results(nn,npp,1,5) = sum(eig0)/n/prod(abs(eig0).^(1/n));
        eig1 = eig(A1);
        results(nn,npp,2,5) = sum(eig1)/n/prod(abs(eig1).^(1/n));
        eign = eig(An);
        results(nn,npp,3,5) = sum(eign)/n/prod(abs(eign).^(1/n));
        eigs = eig(As);
        results(nn,npp,4,5) = sum(eigs)/n/prod(abs(eigs).^(1/n));
        eiga = eig(Aa);
        results(nn,npp,5,5) = sum(eiga)/n/prod(abs(eiga).^(1/n));
        
        
        xinit = zeros(n,1);  
        maxit = 5000;
        fprintf('We go to the solvers: %s \n', solver)

        if strcmp(solver,'lsqr')
            start0 = tic;
            [x0,flag0,relres0,iter0,resvec0,lsvec0] = lsqr(A0,B,tol,maxit,[],[],xinit);
            stop0 = toc(start0);
            disp('Done 0! \n')
            start1 = tic;
            [x1,flag1,relres1,iter1,resvec1,lsvec1] = lsqr(A1,B,tol,maxit,[],[],xinit);
            stop1 = toc(start1);
            disp('Done 1! \n')
            startn = tic;
            [xn,flagn,relresn,itern,resvecn,lsvecn] = lsqr(An,B,tol,maxit,[],[],xinit);
            stopn = toc(startn);
            disp('Done n! \n')
            starts = tic;
            [xs,flags,relress,iters,resvecs,lsvecs] = lsqr(As,B,tol,maxit,[],[],xinit);
            stops = toc(starts);
            disp('Done s! \n')
            starta = tic;
            [xa,flaga,relresa,itera,resveca,lsveca] = lsqr(Aa,B,tol,maxit,[],[],xinit);
            stopa = toc(starta);
            disp('Done a! \n')
        elseif strcmp(solver,'cgs')
            % Here for preconditioned conjugate gradient
            start0 = tic;
            [x0,flag0,relres0,iter0,resvec0] = cgs(A0,B,tol,maxit,[],[],xinit);
            stop0 = toc(start0);
            start1 = tic;
            [x1,flag1,relres1,iter1,resvec1] = cgs(A1,B,tol,maxit,[],[],xinit);
            stop1 = toc(start1);
            startn = tic;
            [xn,flagn,relresn,itern,resvecn] = cgs(An,B,tol,maxit,[],[],xinit);
            stopn = toc(startn);
            starts = tic;
            [xs,flags,relress,iters,resvecs] = cgs(As,B,tol,maxit,[],[],xinit);
            stops = toc(starts);
            starta = tic;
            [xa,flaga,relresa,itera,resveca] = cgs(Aa,B,tol,maxit,[],[],xinit);
            stopa = toc(starta);
        elseif strcmp(solver,'pcg')
            % Here for conjugate gradient
            start0 = tic;   % gamma = 0
            [x0,flag0,relres0,iter0,resvec0] = pcg(A0,B,tol,maxit,[],[],xinit);
            stop0 = toc(start0);
            start1 = tic;   % gamma = 1
            [x1,flag1,relres1,iter1,resvec1] = pcg(A1,B,tol,maxit,[],[],xinit);
            stop1 = toc(start1);
            startn = tic;   % gamma = u^-2
            [xn,flagn,relresn,itern,resvecn] = pcg(An,B,tol,maxit,[],[],xinit);
            stopn = toc(startn);
            starts = tic;   % gamma = gammastar
            [xs,flags,relress,iters,resvecs] = pcg(As,B,tol,maxit,[],[],xinit);
            stops = toc(starts);
            starta = tic;   % gamma = gammastar_approx
            [xa,flaga,relresa,itera,resveca] = pcg(Aa,B,tol,maxit,[],[],xinit);
            stopa = toc(starta);
        end
        
        fprintf('\n\n [n,r,t] = [%i %i %i] problem %i out of %i \n',n,r,t,npp,np)
        fprintf('number of iters resp  %i %i %i %i %i\n',iter0,iter1,itern,iters,itera)
        fprintf('condest #s %g %g %g %g %g\n',condest(A0),condest(A1), condest(An),condest(As),condest(Aa))
        fprintf('residuals.  %g %g %g %g %g\n',relres0,relres1,relresn,relress,relresa)
        fprintf('flags.  %i %i %i %i %i\n',flag0,flag1,flagn,flags,flaga)
        fprintf('ress. norm(Ax-b)  %g %g %g %g\n', ...
	         norm(A0*x0-B)/norm(B),norm(A1*x1-B)/norm(B),norm(As*xs-B)/norm(B),norm(Aa*xa-B)/norm(B))
        %storing results
        % for gamma=0
        results(nn,npp,1,2) = iter0;
        results(nn,npp,1,3) = relres0;
        results(nn,npp,1,4) = stop0;
        results(nn,npp,1,6) = stop0; % total time = solve time + gamma time
        % for gamma=1
        results(nn,npp,2,2) = iter1;
        results(nn,npp,2,3) = relres1;
        results(nn,npp,2,4) = stop1;
        results(nn,npp,2,6) = stop1;
        % for gamma=gamman
        results(nn,npp,3,2) = itern;
        results(nn,npp,3,3) = relresn;
        results(nn,npp,3,4) = stopn; %+ stopcompn;   % only for solve time
        results(nn,npp,3,6) = stopn + stopcompn;
        % for gamma=gammastar
        results(nn,npp,4,2) = iters;
        results(nn,npp,4,3) = relress;
        results(nn,npp,4,4) = stops; %+ min(stopgamma, stopgammaC);
        results(nn,npp,4,6) = stops; %+ min(stopgamma, stopgammaC);
        % for gamma=gammastarA (approx)
        results(nn,npp,5,2) = itera;
        results(nn,npp,5,3) = relresa;
        results(nn,npp,5,4) = stopa; %+ stopgammaA;
        results(nn,npp,5,6) = stopa + stopgammaA;
    end % end for problem fixed dimension
    % means
    mean0 = sum(results(nn,:,1,:))/np;
    mean1 = sum(results(nn,:,2,:))/np;
    meann = sum(results(nn,:,3,:))/np;
    means = sum(results(nn,:,4,:))/np;
    meana = sum(results(nn,:,5,:))/np;
    meantimegamma1 = sum(timegammastar(nn,:,1))/np; % spectral
    meantimegamma2 = sum(timegammastar(nn,:,2))/np; % cholesky
    meantimegamma3 = sum(timegammastar(nn,:,3))/np; % approx
    meantimegamma4 = sum(timegammastar(nn,:,4))/np; % 1/u.^2


    %save for every dimension
    
    resultsn(nn,1,:) = mean0; % 1) condest; 2) iter; 3) relres;
    resultsn(nn,2,:) = mean1; % 4) t.solve; 5) omega cond; 6) total time
    resultsn(nn,3,:) = meann;
    resultsn(nn,4,:) = means;
    resultsn(nn,5,:) = meana;
    if meantimegamma1 < meantimegamma2 % spectral is faster
        timegammastarn(nn,1) = meantimegamma1; % spectral
        timegammastarn(nn,2) = 'S';            % save what decomp is used
    elseif meantimegamma1 > meantimegamma2 % cholesky is faster
        timegammastarn(nn,1) = meantimegamma2; % cholesky
        timegammastarn(nn,2) = 'C';
    else
        timegammastarn(nn,1) = meantimegamma1; % spectral
        timegammastarn(nn,2) = 'B'; % b stands for 'both are the same'
    end
    timegammastarn(nn,3) = meantimegamma3; % approx
    timegammastarn(nn,4) = meantimegamma4; % 1/u.^2
    resultsn(nn,4,6) = resultsn(nn,4,6) + timegammastarn(nn,1);

    %%% Save total time, i.e., solve + gamma time
    % totaltimegammastarn(nn,1) = stops + timegammastarn(nn,1);
    % totaltimegammastarn(nn,2) = stopa + timegammastarn(nn,3);
    % totaltimegammastarn(nn,3) = stopn + timegammastarn(nn,4);
end %end for dimension

fig1 = figure(1);
clf
% semilogy(resvec0,'--',linewidth = 1)
% hold on
% semilogy(resvec1,'-k',linewidth = 1)
% semilogy(resvecs,'-.b',linewidth = 1)
% semilogy(resvecs,'-.r',linewidth = 1)
% legend('\gamma = 0','\gamma = e','\gamma = 1/\|U\|^2', '\gamma = \gamma^{\ast}','location','best')
% title([solver, ' for  A(\gamma) x = c'])
% xlabel('Iterations')
% ylabel('||A(\gamma) x - c||')
% hold off
% figname = ...
%    ['../../../condnumbandQPlatexfiles.d/plotomega', solver, '.png'];
% if options.saveplot == true
%     exportgraphics(fig1, figname)
% end





if options.savetable == true
    tableomega
end




% %%%% Create performance profile
% optionsperf.maxit = maxit;
% optionsperf.tol = tol;
% T = performanceprofilegen(results,timegammastar,optionsperf); %data for performance profile
% perf(T,1,solver)
% 
% figname = ...
%    ['../../../condnumbandQPlatexfiles.d/perfprofomega', solver, '.png'];
% if options.saveplot == true
%     exportgraphics(fig1, figname)
% end


%%%% Create performance profile
optionsperf.maxit = maxit;
optionsperf.tol = tol;
T = performanceprofilegen(results,timegammastar,optionsperf);
perf(T,1,'Time ');
fig1 = gca;
% figname = ...
%    ['../../../condnumbandQPlatexfiles.d/', perfprofnametime, solver, '.png'];
figname = ...
   [perfprofnametime, solver, '.png'];
if options.saveplot == true
    exportgraphics(fig1, figname)
end

%%%% Create performance profile iterations
optionsperf.maxit = maxit;
optionsperf.tol = tol;
Tit = performanceprofilegeniter(results,optionsperf);
perf(Tit,1,'Iterations ')
fig1 = gca;
% figname = ...
%    ['../../../condnumbandQPlatexfiles.d/', perfprofnameiter, solver, '.png'];
figname = ...
   [perfprofnameiter, solver, '.png'];
if options.saveplot == true
    exportgraphics(fig1, figname)
end


% profile report
