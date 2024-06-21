%%% Aug21/23  This file is for generating the table comparing wcond
%%% evaluation for three techniques:
%%%           using exact d  eig(A)  lu(A)   chol(A)
clear
seed=1;
rng(1);
n = [500,1000,2000];  % diff. dimensions for n
savetable = true; %true for generating the table, false otherwise
numbprob = 10;   % # of tests for each dimension/kappa condition number
sfrac = 5;   %fraction of large #s is 1/sfrac
xkappa = 2:9;    % order of the kappa cond numb.
%largens = 1e2;
%smallns = 1e-2;
resultstime = zeros(length(n), length(xkappa),numbprob, 3); % 1 for eig, 2 for chol,
resultsprec = zeros(length(n), length(xkappa), numbprob, 3); % 3 for lu
%% kappa(A) is approx.   largens/smallns
for jj = 1:length(n)
    nn = n(jj);
    for kk = 1:length(xkappa)
        xx = xkappa(kk);
        fprintf('\nusing four methods: exact-d  eig(A)  chol(A)  lu(A)  \n')
        for ii = 1:numbprob  % number of tests
           orderlarg = randi(xx-1);     % order of largens
           ordersmall = orderlarg - xx; % order of smallns
           largens = 10^orderlarg;
           smallns = 10^ordersmall;
           [q,~] = qr(randn(nn));
           d = 10*(rand(nn,1)+.5);  % for eigs
           nt = floor(nn/sfrac);
           d(1:nt) = largens*d(1:nt);   % large nos
           d(nt+1:end) = smallns*d(nt+1:end);   % small nos
           tic
           condAexact = (sum(d)/nn)/(prod(d.^(1/nn))); % exact calculation
           timeexact = toc;
           A = q*diag(d)*q';
           A = (A+A')/2;
           tic
           eigA = eig(A);
           %weigdetn = ((prod(eigA.^(1/n))));  % denominator
           weig = (sum(eigA)/(prod(eigA.^(1/nn))))/nn;
           timeweig = toc;
           tic 
           cholA = chol(A);
           wchol = (trace(A)/prod(diag(cholA).^(2/nn)))/nn;
           timewchol = toc;
           tic 
           [l,u,p] = lu(A);
           wlu  = (trace(A)/prod(abs(diag(u)).^(1/nn)))/nn;
           timewlu = toc;

           fprintf('\n\n')
           fprintf('std cond numb kappa %g \n',cond(A))
           fprintf('times:  [%g %g  %g %g] \n', ...
                          timeexact,timeweig,timewchol,timewlu)
           values =   [condAexact weig wchol wlu]';
           fprintf('omega-values:  [%g %g  %g %g] \n',values)
           diffexact =   values-condAexact;
           fprintf('values-exact:  [%g %g  %g %g] \n',diffexact)
           fprintf('rel. std of four omega values:  %g \n',std(values)/mean(values))
           
           %%% saving the resurts
           resultstime(jj,kk,ii,1) = timeweig;
           resultstime(jj,kk,ii,2) = timewchol;
           resultstime(jj,kk,ii,3) = timewlu;
           resultsprec(jj,kk,ii,1) = abs(diffexact(2));
           resultsprec(jj,kk,ii,2) = abs(diffexact(3));
           resultsprec(jj,kk,ii,3) = abs(diffexact(4));
        end  % of for 
    end  % of for xkappa
end  % of for n
meanresulttimes = sum(resultstime,3)/length(numbprob);
meanresultprec = sum(resultsprec,3)/length(numbprob);

if savetable == true
    tableomegatimes
    tableomegaprec
end

