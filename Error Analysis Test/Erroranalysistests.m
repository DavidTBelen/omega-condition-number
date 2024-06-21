%%% file Revisionerroranalysistests.m for condestimates/kappa/omega
%%% revising Nov12 
%https://www.mathworks.com/matlabcentral/fileexchange/97267-numerical-differentiation-toolbox %added
%% doc_NDT     or    doc_NDT name  e.g.  doc_NDT cderivative
%%% initial data n
clear
%addpath(genpath('numdiff.d'))
seed = 100;
rng(seed);
%saverng = rng('shuffle');
%save('savefile','saverng');
n = 200;
%%% comment out the distributions NOT used
kind = 1; % uniform distribution for eigs
%kind = 2; % normal shifted distribution for eigs
if kind == 1
        fprintf('using unif. distrib eigs: deigs = 10*(rand(n,1)+.1)\n');
elseif kind == 2
        fprintf('using unif. distrib eigs: deigs = 10*(rand(n,1)+.1)\n');
end
density = .1;   % for sparse A pos def
ntestsbbar = 1000;   % bbar and db input to black box/dir. deriv.
perteps = 1e-6;    % accuracy for derivatives
kappas = (1:200)';    % target cond(A) lambdamax/lambdamin
numbAs = length(kappas);
omega = zeros(numbAs,1);
condepsm = zeros(numbAs,1);  

%%%%%%%%%%%%%%%%% for loop for different As
In = eye(n);
ei = In(:,randi([1 n],1)); % choose a random index for partial deriv.
for jj = 1:numbAs  % different matrices A_jj
   %kappa = kappas(jj);  % desired cond(A)
   [VA,~] = qr(randn(n));
   if kind == 1  % unif distrib eigs
           deigs = (rand(n,1)+.1);% uniform distrbn
           A = VA*diag(deigs)*VA';
           A = (A+A')/2;
	   if min(eig(A)) < 1e-6
		   fprintf('not pos def\n')
		   keyboard
	   end
   elseif kind == 2  % use normal distrn eigs
           deigs = randn(n,1);  % normal distribution
	   if min(deigs) < 0
	           deigs = deigs - min(deigs) + .1;
	   end
           A = VA*diag(deigs)*VA';
           A = (A+A')/2;
	   if min(eig(A)) < 1e-6
		   fprintf('not pos def\n')
		   keyboard
	   end
   end
   kappas(jj) = cond(A);
   cholA = chol(A);
   omega(jj) = (trace(A)/prod(diag(cholA).^(2/n)))/n;
   condeps = zeros(ntestsbbar ,1);  
   %% given  A fixed; input bbar, direction db  
   %% evaluate dir derivative of  f_A'(bbar,db)
   for ii = 1:ntestsbbar % number of tests for bbar and db
      %xbar = VA(:,round(n/2)-2:round(n/2)+2)*randn(5,1);
      xbar = randn(n,1);
      bbar = A*xbar;  % exact solution
      db = randn(n,1);    % perteps = 1e-6;    ???
      db = perteps*db/norm(db);
      %rhob = perteps*norm(db)/norm(bbar);
      %f = @(bbar,db,A,t)(norm(A\(bbar+t*db)));
      %f = @(t)(norm(A\(bbar+t*db)));
      %f = @(t)(abs(ei'*(A\(bbar+t*db))));
      %condeps(ii) = abs(cderivative(f,0));
      xdx = A\db; %(bbar+perteps*db);
      condeps(ii) = (norm(xdx)*norm(bbar))/(norm(xbar)*norm(db));  %(norm(xdx-xbar)/norm(xbar))/rhob;
      if condeps(ii) > kappas(jj)
   	   fprintf('condeps > kappa??\n')
   	   keyboard
      end

   end  % of for random bbar, db
   condepsm(jj) = mean(condeps);
end  % of for   A pos def,b


%%%%%%%%%%%%%%%%linear models
   m = numbAs;
   X = [ones(m,1) condepsm];
   betak = X\kappas;
   ykappa = betak(1) + betak(2)*condepsm;
   X = [ones(m,1) condepsm];
   betao = X\omega;
   yomega = betao(1) + betao(2)*condepsm;
   h1 = figure(1);
   clf
   box on
   hold on
   scatter(condepsm, kappas,'LineWidth',1.4)
   plot(condepsm,ykappa,'LineWidth',1.4)
%    plot(1:m,yomega)
%    plot(1:m,condepsm)
   legend('(cond(A_i),\kappa(A_i))','Regression line','location','southeast')
   if kind == 1
    title('LRM cond and \kappa(A) for uniform distributed eigenvalues')
   elseif kind == 2
    title('LRM cond and \kappa(A) for normally distributed eigenvalues')
   end
   xlabel('Mean of  condition number estimates')
   ylabel('\kappa(A)')
   xlim([min(condepsm),max(condepsm)])
   set(gca,'ClippingStyle','rectangle');
   %set(gca,'LooseInset',get(gca,'TightInset')); 
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
   if kind == 1
    print(gcf,'..\Latexrevisingasofnov9COAPSubmission\lrm_kappa_uniform.pdf','-dpdf','-r300','-bestfit');
   elseif kind == 2
    print(h1,'..\Latexrevisingasofnov9COAPSubmission\lrm_kappa_normal.pdf','-dpdf','-r300','-bestfit');
   end
   hold off

   h2 = figure(2);
   clf
   box on
   hold on
   scatter(condepsm, omega,'LineWidth',1.4)
   plot(condepsm,yomega,'LineWidth',1.4)
   legend('(cond(A_{i}),\omega(A_{i}))','Regression line','location','southeast')
   if kind == 1
    title('LRM cond and \omega(A) for uniform distributed eigenvalues')
   elseif kind == 2
    title('LRM  cond and \omega(A) for normally distributed eigenvalues')
   end
   xlabel('Mean of condition number estimates')
   ylabel('\omega(A)')
   xlim([min(condepsm),max(condepsm)])
    %set(gca,'LooseInset',get(gca,'TightInset')); 
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
   if kind == 1
    print(gcf,'..\Latexrevisingasofnov9COAPSubmission\lrm_omega_uniform.pdf','-dpdf','-r300','-bestfit');
   elseif kind == 2
    print(h2,'..\Latexrevisingasofnov9COAPSubmission\lrm_omega_normal.pdf','-dpdf','-r300','-bestfit');
   end
   hold off



%    figure(2)
%    clf
%    hold on
%    plot(1:m,abs(ykappa-condepsm))
%    plot(1:m,abs(yomega-condepsm))
%    legend('abs(ykappa-condepsm)','abs(yomega-condepsm)', ...
% 	           'location','best')
%    title('absolute values of the diffences for the linear models')
%    hold off
%    drawnow
%    
corrkappa =  corrcoef(kappas,condepsm);
corromega = corrcoef(omega,condepsm);
fprintf('correlation coeff kappa to condepsm, omega to condepsm %g %g\n', ...
	      corrkappa(1,2),corromega(1,2))


%%%%%%%%%%%%%%%%quadratic models
%    m = numbAs;
%    X = [ones(m,1) kappas kappas.^2];
%    beta = X\condepsm;
%    ykappa = beta(1) + beta(2)*kappas + beta(3)*kappas.^2;
%    X = [ones(m,1) omega omega.^2];
%    beta = X\condepsm;
%    yomega = beta(1) + beta(2)*omega + beta(3)*omega.^2;
%    figure(3)
%    clf
%    hold on
%    plot(1:m,ykappa)
%    plot(1:m,yomega)
%    plot(1:m,condepsm)
%    legend('kappa','omega','condepsm','location','best')
%    title('two quadratic regression models: omega/kappa')
%    hold off
% 
% 
% 
%    figure(4)
%    clf
%    hold on
%    plot(1:m,abs(ykappa-condepsm))
%    plot(1:m,abs(yomega-condepsm))
%    legend('abs(ykappa-condepsm)','abs(yomega-condepsm)', ...
% 	           'location','best')
%    title('absolute values of the diffences for the quadratic models')
%    hold off
%    drawnow
   
