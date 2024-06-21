%%% This file generates the table and performance profile in Section 4.2.
%%% 
clear
seed = 1;
rng(seed);
profile clear
profile off

options.savetable = true; %true for save table
options.saveplot = true; %true for sace plot


sizesn = [100,200,500,1000,2000];  %dimensions for tests
resultsn = zeros(length(sizesn),3,5);
optionsplot.sizesn = sizesn;

solver = 'pcg'; % 'lsqr' or 'cgs'
optionsplot.solver =  solver;

timegammastarn = zeros(length(sizesn),2);

%%% table name and performance profile figure name
tablename = 'tableomega';   % name used for latex file as of aug15/23
perfprofnameiter = 'perfprofiter';  % name used for latex file as of aug15/23
perfprofnametime = 'perfprofomega';  % name used for latex file as of aug15/23
%  ['../../../condnumbandQPlatexfiles.d/perfprofiter', solver, '.png'];


np = 10; %number of problems
%np = 5; %number of problems
%np = 2; %number of problems

setting = zeros(length(sizesn),np,4); %for storing n,r and t
results = zeros(length(sizesn),np,4,5);  %for storing iterations, cond, residuals and time
timegammastar = zeros(length(sizesn),np,2); %for storing the time for computing gammastar
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

        Aeps = A + (max(0,-min(eig(A)))+epsilon)*speye(n);
        A0 = Aeps;
        A1 = Aeps + U*U';
        startcompn = tic;
        nun = full((norms(U).^2)');
        gamman = 1./nun;
        gamman(gamman>1) = 1;
        stopcompn = toc(startcompn);
        An = Aeps + U*diag(gamman)*U';
        An = (An+An')/2;
        As = Aeps + U*diag(gammastar)*U';
        As = (As+As')/2;
        
        % storing the condest condition number
        results(nn,npp,1,1) = condest(A0);
        results(nn,npp,2,1) = condest(A1);
        results(nn,npp,3,1) = condest(An);
        results(nn,npp,4,1) = condest(As);
        % storing the omega condition number
        eig0 = eig(A0);
        results(nn,npp,1,5) = sum(eig0)/n/prod(abs(eig0).^(1/n));
        eig1 = eig(A1);
        results(nn,npp,2,5) = sum(eig1)/n/prod(abs(eig1).^(1/n));
        eign = eig(An);
        results(nn,npp,3,5) = sum(eign)/n/prod(abs(eign).^(1/n));
        eigs = eig(As);
        results(nn,npp,4,5) = sum(eigs)/n/prod(abs(eigs).^(1/n));
        
        
        xinit = zeros(n,1);  
        maxit = 5000;
        disp('We go to the solvers: \n')

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
        elseif strcmp(solver,'pcg')
            % Here for conjugate gradient
            start0 = tic;
            [x0,flag0,relres0,iter0,resvec0] = pcg(A0,B,tol,maxit,[],[],xinit);
            stop0 = toc(start0);
            start1 = tic;
            [x1,flag1,relres1,iter1,resvec1] = pcg(A1,B,tol,maxit,[],[],xinit);
            stop1 = toc(start1);
            startn = tic;
            [xn,flagn,relresn,itern,resvecn] = pcg(An,B,tol,maxit,[],[],xinit);
            stopn = toc(startn);
            starts = tic;
            [xs,flags,relress,iters,resvecs] = pcg(As,B,tol,maxit,[],[],xinit);
            stops = toc(starts);
        end
        
        fprintf('\n\n [n,r,t] = [%i %i %i]\n',n,r,t)
        fprintf('number of iters resp  %i %i %i %i\n',iter0,iter1,itern,iters)
        fprintf('condest #s %g %g %g %g\n',condest(A0),condest(A1), condest(An),condest(As))
        fprintf('residuals.  %g %g %g %g\n',relres0,relres1,relresn,relress)
        fprintf('flags.  %i %i %i %i\n',flag0,flag1,flagn, flags)
        fprintf('ress. norm(Ax-b)  %g %g %g\n', ...
	         norm(A0*x0-B)/norm(B),norm(A1*x1-B)/norm(B),norm(As*xs-B)/norm(B))
        %storing results
        % for gamma=0
        results(nn,npp,1,2) = iter0;
        results(nn,npp,1,3) = relres0;
        results(nn,npp,1,4) = stop0;
        % for gamma=1
        results(nn,npp,2,2) = iter1;
        results(nn,npp,2,3) = relres1;
        results(nn,npp,2,4) = stop1;
        % for gamma=gammastar
        results(nn,npp,3,2) = itern;
        results(nn,npp,3,3) = relresn;
        results(nn,npp,3,4) = stopn + min(stopcompn,stopgammaC);
        % for gamma=gammastar
        results(nn,npp,4,2) = iters;
        results(nn,npp,4,3) = relress;
        results(nn,npp,4,4) = stops;
    end % end for problem fixed dimension
    % means
    mean0 = sum(results(nn,:,1,:))/np;
    mean1 = sum(results(nn,:,2,:))/np;
    meann = sum(results(nn,:,3,:))/np;
    means = sum(results(nn,:,4,:))/np;
    meantimegamma1 = sum(timegammastar(nn,:,1))/np;
    meantimegamma2 = sum(timegammastar(nn,:,2))/np;


    %save for every dimension
    
    resultsn(nn,1,:) = mean0;
    resultsn(nn,2,:) = mean1;
    resultsn(nn,3,:) = meann;
    resultsn(nn,4,:) = means;
    timegammastarn(nn,1) = meantimegamma1;
    timegammastarn(nn,2) = meantimegamma2;

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
perf(T,1,'Time ')

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

% figname = ...
%    ['../../../condnumbandQPlatexfiles.d/', perfprofnameiter, solver, '.png'];
figname = ...
   [perfprofnameiter, solver, '.png'];
if options.saveplot == true
    exportgraphics(fig1, figname)
end


profile report
