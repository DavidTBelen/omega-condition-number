%%%minimize use Dindelbach to solve fractional progr for opt cond. #
%% as of dec18/24  add comparison  to new opt dbar for Ainv opt diag precond
%% as of oct30/24 this file draws figure on improved clustering of
%% eigenvalues in the paper in section 2; currently figure 2.1.
%%% i.e. figure 3 below is figure 2.1 in the paper.
%%% we are comparing opt diagonal scalings    D W D on W
%file is   optprecondomegakappa.m
clear
%seed = randi(20,1);
%seed = 11;
%seed = 1430034934;  % nice pic???
%fprintf('seed is %i\n',seed)
%rng(seed);
seed = 52124;  % seed 52124 for paper as of oct30/24 figure 3 is
                          %% figure 2.1 in the paper about clustering
rng(seed)
%saverng = rng('shuffle');
%%%%%%%Initialization with random data
%%%% first choose which condition number
kappatf = true;  % for kappa cond #; else false means kappa not in figure
%if kappatf
%	omegatf = false;   % ???need??? not used???
%else
%	omegatf = true;
%end
randS = false;
if randS
	n = 80;
	density = .5;
	rc = .01;
	%S = 10*sprandsym(n,density,rc,1);  % role of V in prop. in paper
	%W = S^2;
	W = 10*sprandsym(n,density,rc,1);
	S = sqrtm(full(W));
	S(abs(S)<1e-10) = 0;
	S = (S + S')/2;
	figure(1)
	clf
	spy(S)
	title('Matrix S')
	fprintf('min eig S %g \n',min(eig(S)))
	S = sparse(S);
else
	n = 30;
	density = .8;
	rc = .1;
	%S1 = sprandsym(n,density,rc)+diag(10*rand(n,1));
	S1 = sprandsym(n,density,rc);
	S2 = 300*(sprandsym(n,density,rc)+ max(eig(S1))*speye(n));
	S12 = .1*sprandn(n,n,.7);
	S = [S1 S12; S12' S2];
	n = length(S);
	W = S^2;
	figure(1)
	clf
	spy(W)
	title('Matrix W')
end

%dhat = (1/n)*ones(n,1);  %starting vector
%%% we are finding D^2 for opt scaling D using Dinkelbach
dhat = rand(n,1);  %starting vector  for Dinkelbach
dhat = dhat/norm(dhat);
dhat = norm(W,'fro')*dhat/2;
omega = @(Y)( (trace(Y)/n)/(det_rootn(Y))  );
YRinv = @(Y)( inv(chol(Y)) );  % Y = chol(Y)'*chol(Y); Yinv = YRinv*YRinv'
Winvsq = @(Y)( (YRinv(Y)*YRinv(Y)')^2 );       %  (WRinv*WRinv')^2;
omegaAinv =  @(Y)( sqrt(omega(Winvsq(Y))) );
%%%% D=diag(d) ;  S D S == sqrt D W sqrt D  ; eigs and kappa,omega unchanged 
Sd = @(d)( (S*diag(d)*S) );   % Wdiag(d) = Sdiag(d)S for tr/det/eig fns
Wd = @(d)( (diag(d)*W*diag(d)) );   % for dnew opt scaling
fomega = @(d)( (trace(Wd(d)))/n );  %  for DWD
gomega = @(d)(  det_rootn(Wd(d)) );  % can make MORE efficient????
Fomega = @(d)(  fomega(d)/gomega(d) ); % for DWD
fkappa = @(d)( ( lambda_max(Sd(d))  )  );  % for S DD S; d == DD
gkappa =  @(d)( (lambda_min(Sd(d)))  );  % for S DD S; d == DD
Fkappa = @(d)(  fkappa(d)/gkappa(d) );
qhatomega = Fomega(dhat);  % initial data
qhatkappa = Fkappa(dhat);  % initial data
fkappaW = @(W)( ( lambda_max(W)  )  );
gkappaW =  @(W)( (lambda_min(W))  );
FkappaW = @(W)(  fkappaW(W)/gkappaW(W) );

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

%%%% use Dinkelbach to find D*D
tic
while Fq  < -reltoler && iter <= maxiter

	cvx_begin  quiet
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
%%%% ds  SQUARED optimal diagonal scaling for kappa
dsfinal = ds;
%%% three opt precond folllow
ds = sqrt(ds);  % opt scaling is diag(ds)*W*diag(ds) for kappa
%%%%% opt. omega scale:
domegaopt = 1./sqrt(diag(W));  % theoretical opt precond. for omega
%%%% add new preconditioner
dnewAinv = testNewtonfornewdbar(W);
dbarnewomegaopt = 1./dnewAinv;   % for new Ainv^2 omega
%%%%%%%%%%%%%
%%% why next line??? why scale 'nicely'????
% domegaopt = (ds(1)/domegaopt(1))*domegaopt;  % scale 'nicely'
fprintf('End of Dinkelbach: time = %g;  final Fq = %g\n',timewhile,Fq)
omegaW = omega(W);
Womega = Wd(domegaopt);   % opt omega scaling
Womega = (Womega+Womega')/2;
omegaWomega = omega(Womega);
WomegaAinv = Wd(dbarnewomegaopt);   % opt omega Ainvsqrd scaling
WomegaAinv = (WomegaAinv+WomegaAinv')/2;
omegaWAinv = omega(WomegaAinv); 
Wkappa = Wd(ds);   % Dkappa W Dkappa
Wkappa = (Wkappa+Wkappa')/2;
omegaWkappa= omega(Wkappa);

eigsW = eig(W);
eigsWomega = eig(Womega);
eigsWomegaAinv = eig(WomegaAinv);
eigsinvWomegaAinv = 1./eigsWomegaAinv;
eigsWkappa = eig(Wkappa);


if kappatf
	fprintf('KAPPA cond vals: ds, domegaopt, origW :  %g %g %g \n', ...
		CondF(ds),CondF(domegaopt),CondF(ones(n,1)))
	if norm(ds-domegaopt) > 1e-9
		fprintf('ERROR for kappa code: ds - domegaopt?? %g\n', ...
			norm(ds-domegaopt))
	end
	figure(3)
	clf
	semilogy(eigsW)
	hold on
	semilogy(eigsWomega)
	semilogy(eigsWomegaAinv)
	semilogy(eigsWkappa)
	semilogy(sort(eigsinvWomegaAinv,'ascend'));
	legend('W','omega  opt','omega A inverse opt','kappa opt', ...
		'reciprocal omega A inverse opt ', ...
		'location','best')
	title('eigenvalues distributions for W and preconditioned W','Fontsize',15)
	xlabel('ith eigenvalue','Fontsize',15)
	ylabel('magnitude of eigenvalue','Fontsize',15)
	hold off
	%% find the std for middle eigs? i.e.  best 3 sets
	stdW = std(eigsW);  % orig eigs W
	stdWomega = std(eigsWomega);  % omega opt eigs W
	stdWomegaAinv = std(eigsWomegaAinv);  % omega opt eigs W^-2 new d
	stdWkappa = std(eigsWkappa);  % kappa opt eigs W
	%%% in future find the best three set clustering!!!?????
	%%% for now use 20 points at each end
	%%This is for fig 2.1 comment as of Oct30/24
	fprintf('std eigs: W, Womega, Wdnew, Wkappa is: %g %g %g %g \n', ...
	        stdW,stdWomega,stdWomegaAinv,stdWkappa)
	fprintf('std 1:end-45 eigs for W, Womega, Wnew, Wkappa is: %g %g %g %g \n', ...
	  std(eigsW(1:end-45)),std(eigsWomega(1:end-45)), ...
	          std(eigsWomegaAinv(1:end-45)),std(eigsWkappa(1:end-45)))
	fprintf('omega of: W,Womega,Womeganew,Wkappa: %g %g %g %g \n', ...
               omegaW,omegaWomega,omegaWAinv,omegaWkappa)
	fprintf('kappa of: W,Womega,Womeganew,Wkappa: %g %g %g %g \n', ...
		FkappaW(W),FkappaW(Womega),FkappaW(WomegaAinv),FkappaW(Wkappa))
	fprintf('omegaAinv of: W,Womega,Womeganew,Wkappa: %g %g %g %g \n', ...
          omegaAinv(W),omegaAinv(Womega), ...
	           omegaAinv(WomegaAinv),omegaAinv(Wkappa))

	   %%% mean of eigs of DWD is 1

	meanseigs(1) =  mean(eigsW);
	meanseigs(2) = mean(eigsWomega);
	meanseigs(3) = mean(eigsWomegaAinv);
	meanseigs(4) = mean(eigsWkappa);
	eigsone = zeros(4,5);
	iis = -1:.1:3;
	for ii = 1:length(iis)
		eigsone(1,ii) = sum(abs(meanseigs(1)-eigsW)<10^iis(ii));
		eigsone(2,ii) = sum(abs(meanseigs(2)-eigsWomega)<10^iis(ii));
		eigsone(3,ii) =  ...
			 sum(abs(meanseigs(3)-eigsWomegaAinv)<10^iis(ii));
		eigsone(4,ii) = ...
			sum(abs(meanseigs(4)-eigsWkappa)<10^iis(ii));
	end
	figure(4)
	clf
	plot(iis(1:end),eigsone(1,:),'x-')
	hold on
	plot(iis(1:end),eigsone(2,:),'o-')
	plot(iis(1:end),eigsone(3,:),'d-')
	plot(iis(1:end),eigsone(4,:),'s-')
	legend('eigsW','eigsWomega','eigsomegaAinv','eigskappa', ...
		'location','best')
	title('number of eigs within decimals from the mean','Fontsize',15)
	xlabel('order of width of band about the mean','Fontsize',15)
	ylabel('number of eigs within band','Fontsize',15)
	hold off

	        

else
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End Dinkelbach algor.



