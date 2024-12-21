%function comparecondnumbofcondnumb
%%%% Henry modifying Nov3/24
% use function to help with profile ????
%%% compare condition number calc times and values
%%  note that density of A greatly changes the comparison as
%% omega does much better with sparse problems
%% see the  plot for what this file is doing
clear
profilerep = false;
if profilerep
	profile clear
	profile on
end
density0 = 1e-4;
rcs = .1.^(linspace(1,3,10));
%rcs = .1.^[2:6];
noisemag = 1e-8;   % magnitude to perturb eigs
n = 1000;
midpteigs = true; %%    for mid points multiple eigs
if midpteigs
	n = 2*floor(n/2)+1;   % ensure it is odd
	midn = floor(n/2)+1;   % midpoint
	fprintf('Prob: increasing multiplicity of midpoint eigenvalue; ')
	fprintf('  n = %i\n')
end
density = density0/n;
fprintf('magnitude of noise=%g;  density0/n=%g\n',noisemag,density)

omegaA = @(A)( (trace(A)/n)/det_rootn(A) );

iterrc = 0;
for rc = rcs
iterrc = iterrc +1;
A = sprandsym(n,density,rc,1);
[Q,D] = eig(full(A));
d = diag(D);
dorig = d;
ekappa = d(end)/d(1);  % original exact kappa 
fprintf('with rc = %g:  original kappa and omega are %e %e  \n', ...
	 ekappa,omegaA(A),rc)
dnoise = d;  % recover A using dnoise eigenvalues
iter = 0;
%while iter < floor(n/20);
while iter < 10
	iter = iter + 1;
	A = Q*diag(dnoise)*Q';  % recover perturbed A after eig perturbation
	eeigsA = eig(full(A));
	akappa = eeigsA(end)/eeigsA(1);  % approx kappa
	Kdiff(iter) = abs(ekappa-akappa)/(ekappa+1);
	eomega = mean(d)/geo_mean(d);  % exact omega
        aomega = omegaA(A);  % approx omega
	Odiff(iter) = abs(eomega-aomega)/(eomega+1);
	akappas(iter) = akappa;
	aomegas(iter) = aomega;
	ekappas(iter) = ekappa;  % does not change
	eomegas(iter) = eomega;
	if midpteigs %%    for mid points multiple eigs; pert all eigs
		%d(midn+iter) = d(midn);
		%d(midn-iter) = d(midn);
		noise = randn(n,1);
		noise = noisemag*noise/norm(noise);
		dnoise = d;
		dnoise = dnoise + noise;
	else %%    for end points multiple eigs
		d(1+iter) = d(1);
		d(end-iter) = d(n);
		noise = randn(n-2*iter,1);
		noise = noisemag*noise/norm(noise);
		dnoise = d;
		dnoise(1+iter:end-iter) = ...
			dnoise(1+iter:end-iter)+noise;  %perturb
	end


	%Kdiffr(iterrc) = Kdiff(iter);
	%Odiffr(iterrc) = Odiff(iter);
end      % while iter < 
	Kdiffr(iterrc) = mean(Kdiff);
	Odiffr(iterrc) = mean(Odiff);
	figure(1)
	clf
	plot(rcs(1:iterrc),max(1e-16,Kdiffr(1:iterrc)),'-x')
	hold on
	plot(rcs(1:iterrc),Odiffr(1:iterrc),'-o')
	texttitle = ['cond.# relerr |appr.-exct| vs  rc '];
	title(texttitle)
	if midpteigs
		xlabel('increasing reciprocal condition number')
		ylabel('plot of relerr |appr.-exct|') 
	else
		xlabel('increasing multiplicity of max-min eigenvalues')
	end
	legend('rel. |exact-approx| \kappa','rel. |exact-approx| \omega', ...
		         'location','best')
	hold off
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
    %print(gcf,'CMPplotcondsvsrc.pdf','-dpdf','-r300','-bestfit');
	drawnow
	figure(2)
	clf
	semilogy(rcs(1:iterrc),max(1e-16,Kdiffr(1:iterrc)),'x')
	hold on
	semilogy(rcs(1:iterrc),Odiffr(1:iterrc),'o')
	texttitle = ['cond.# relerr |appr.-exct| vs  rc '];
	title(texttitle)
	if midpteigs
		xlabel('increasing reciprocal condition number')
		ylabel('semilogy of relerr |appr.-exct|') 
	else
		xlabel('increasing multiplicity of max-min eigenvalues')
	end
	legend('rel. |exact-approx| \kappa','rel. |exact-approx| \omega', ...
		         'location','best')
	hold off
	set(gca,'units','centimeters')
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
    %print(gcf,'CMPsemilogycondsvsrc.pdf','-dpdf','-r300','-bestfit');
end      % for rcs
if profilerep
	profile report
end


