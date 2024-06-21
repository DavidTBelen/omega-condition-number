%%%minimize use Dindelbach to solve fractional progr for opt cond. #
%file is   optprecondomegakappa.m
clear
%seed = randi(20,1);
%seed = 11;
%fprintf('seed is %i\n',seed)
%rng(seed);
seed = 52124;  % for paper
rng(seed)
%%%%%%%Initialization with random data
%%%% first choose which condition number
kappatf = true;  % for kappa cond #; else false means omega chosen
if kappatf
	omegatf = false;
else
	omegatf = true;
end
randS = false;
if randS
n = 120;
%W = randn(n);
%W = W'*W;   % r by r pos def
%W = (W'+W)/2;   % r by r pos def
density = .05;
rc = .01;
S = sprandsym(n,density,rc);
W = S^2;
figure(1)
clf
spy(W)
title('Matrix W')
else
n = 80;
density = .05;
rc = .01;
%S1 = sprandsym(n,density,rc)+diag(10*rand(n,1));
S1 = sprandsym(n,density,rc);
S2 = 10*(sprandsym(n,density,rc)+ max(eig(S1))*speye(n));
S12 = .001*sprandn(n,n,.001);
S = [S1 S12; S12' S2];
[n,n] = size(S);
W = S^2;
figure(1)
clf
spy(W)
title('Matrix W')
end

%dhat = (1/n)*ones(n,1);  %starting vector
dhat = rand(n,1);  %starting vector
dhat = dhat/norm(dhat);
dhat = norm(W,'fro')*dhat/2;
%S = sqrtm(W);
omega = @(W)( (trace(W)/n)/(det_rootn(W))  );
Wd = @(d)( (S*diag(d)*S) );   % Wdiag(d) = Sdiag(d)S for tr/det/eig fns
fomega = @(d)( (trace(Wd(d)))/n );
gomega = @(d)(  det_rootn(diag(d))*(det(W^(1/n)))  );
Fomega = @(d)(  fomega(d)/gomega(d) );
%gomega = @(d)(  det_rootn(diag(d))  );   % det W constant
fkappa = @(d)( ( lambda_max(Wd(d))  )  );
gkappa =  @(d)( (lambda_min(Wd(d)))  );
Fkappa = @(d)(  fkappa(d)/gkappa(d) );
qhatomega = Fomega(dhat);  % initial data
qhatkappa = Fkappa(dhat);  % initial data

if kappatf
	f = fkappa;
	g = gkappa;
	CondF = Fkappa;
	qhat = qhatkappa;
else
	f = fomega;
	g = gomega;
	CondF = Fomega;
	qhat = qhatomega;
end


%%%%%%%%find the optimal q using Dinkelbach
dsnew = dhat;
Fq = -inf;    % minimization so F is increasing and is negative
toler = 1e-12;
maxiter = 10;
iter = 1;
Fqs = [];
qs = [];
reltoler =  toler*(norm(W,'fro'));
qs(iter) = CondF(dsnew) + 10;   % ensure f(dsnew) - qs*gsnew < 0 as f>0,g>0
M = 10;  % to start negative
Fqs(iter) =  f(dsnew) - qs(iter)*g(dsnew) - M;

tic
while Fq  < -reltoler && iter <= maxiter

	cvx_begin
	cvx_precision best
	variable ds(n,1) nonnegative
	dual variable lam
	minimize( f(ds) - qs(iter)*g(ds));
	subject to
		lam : norm(ds) <= norm(W,'fro');
	cvx_end


	dsnew = ds;
	Fq = cvx_optval;
	figure(2)
	clf
	%% use -Fqs for semilogy that needs Pos els
	plot(qs(1:iter),Fqs(1:iter),'-x')
	%semilogy(qs(1:iter),-min(0,Fqs(1:iter)),'-x')
	%legend('q value','Fq value','location','best')
	xlabel('q values')
	ylabel('-F(q) values')
	title(['q vs -Fq in ',num2str(iter), ...
		' iters; Dinkelbach alg. for Xsp'])
	drawnow
	iter = iter + 1;
	Fqs(iter) = Fq;
	qs(iter) = CondF(dsnew);
end
timewhile = toc;
%%%%% opt. omega scale:
domegaopt = 1./diag(W);  % theoretical opt precond.
domegaopt = (ds(1)/domegaopt(1))*domegaopt;  % scale 'nicely'
fprintf('End of Dinkelbach: time = %g;  final Fq = %g\n',timewhile,Fq)
figure(3)
eigsW = sort(eig(W));
eigsWd = sort(eig(Wd(domegaopt)));
if omegatf
	fprintf('OMEGA  cond vals: ds,  domegaopt, origW :  %g %g %g \n', ...
		CondF(ds),CondF(domegaopt),CondF(ones(n,1)))
	if norm(ds-domegaopt) > 1e-9
		fprintf('ERROR for omega code: dso - domegaopt?? %g\n', ...
			norm(ds-domegaopt))
	end
	figure(3)
	clf
	semilogy(eigsW)
	hold on
	semilogy(eigsWd)
	legend('eigs W','eigs d omegaopt','location','best')
	title('sorted eigenvalues for W and for omega opt scaled W')
	hold off
end
if kappatf
	fprintf('KAPPA cond vals: ds, domegaopt, origW :  %g %g %g \n', ...
		CondF(ds),CondF(domegaopt),CondF(ones(n,1)))
	if norm(ds-domegaopt) > 1e-9
		fprintf('ERROR for kappa code: ds - domegaopt?? %g\n', ...
			norm(ds-domegaopt))
	end
	eigsds = sort(eig(Wd(ds)));
	figure(3)
	clf
	semilogy(eigsW)
	hold on
	semilogy(eigsWd)
	semilogy(eigsds)
	legend('eigs W','eigs d omegaopt','eigs d kappaopt','location','best')
	title('sorted eigenvalues for: W, omega opt; kappa opt scaled W')
	hold off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End Dinkelbach algor.



