clear all

n = 4000;

r = 5;
Q = zeros(n);
for ii=1:n/r
    tempA = rand(r);
    tempA = (tempA+tempA)'/2;
    
    [tempor,~] = qr(tempA);
    Q(r*(ii-1)+1:r*ii,r*(ii-1)+1:r*ii) = tempor;

 
end

permQ = randperm(n);
Q = Q(permQ,permQ);

Q = sparse(Q);


% problemm = 'nos4';
% load(problemm)
% W = Problem.A;
% [Q,~] = eig(full(W));
% n = length(Q);


lam1 = 4.042307e+11;
lamn = 7.327514e-06;
kappa = lam1/lamn;




% Matrix 1: Uniformly distributed eigenvalues:
d1_0 = (lam1-lamn)*rand(n-2,1)-lamn;
meaneigs1 = sum(d1_0)/(n-2);
d1 = [lamn;d1_0;lam1];
if min(d1)<lamn || max(d1)>lam1
    print('Error with eigenvalues!!')
    keyboard
end
m1 = mean(d1);
std1 = std(d1);
omega1 = sum(d1)/n/(prod(d1.^(1/n)));
M1 = Q*diag(d1)*Q';
M1(abs(M1)< 1e-12) = 0;  % getrid of small roundoff errors
M1 = sparse((M1+M1')/2);


% Matrix 2: Normally distributed eigenvalues:
d2_0 = (lam1-lamn)/10*randn(n-2,1)+(lam1+lamn)/2;
d2_0 = d2_0 + (meaneigs1 -mean(d2_0));
%d2_0 = d2_0 + m1-(lam1+lamn)/n;
d2 = [lamn;d2_0;lam1];
if min(d2)<lamn || max(d2)>lam1
    print('Error with eigenvalues!!')
    keyboard
end
m2 = mean(d2);
std2 = std(d2);
omega2 = sum(d2)/n/(prod(d2.^(1/n)));
M2 = Q*diag(d2)*Q';
M2(abs(M2)< 1e-12) = 0;  % getrid of small roundoff errors
M2 = sparse((M2+M2')/2);

% Matrix 3: Normally distributed with smaller std:
d3_0 = (lam1-lamn)/20*randn(n-2,1)+(lam1+lamn)/2;
d3_0 = d3_0 + (meaneigs1 -mean(d3_0));

%d2_0 = d2_0 + m1-(lam1+lamn)/n;
d3 = [lamn;d3_0;lam1];
if min(d3)<lamn || max(d3)>lam1
    disp('Error with eigenvalues!!')
    keyboard
end
m3 = mean(d3);
std3 = std(d3);
omega3 = sum(d3)/n/(prod(d3.^(1/n)));
M3 = Q*diag(d3)*Q';
M3(abs(M3)< 1e-12) = 0;  % getrid of small roundoff errors
M3 = sparse((M3+M3')/2);


% Matrix 4: Normally distributed with even smaller std:
d4_0 = (lam1-lamn)/50*randn(n-2,1)+(lam1+lamn)/2;
d4_0 = d4_0 + (meaneigs1 -mean(d4_0));

%d2_0 = d2_0 + m1-(lam1+lamn)/n;
d4 = [lamn;d4_0;lam1];
if min(d4)<lamn || max(d4)>lam1
    disp('Error with eigenvalues!!')
    keyboard
end
m4 = mean(d4);
std4 = std(d4);
omega4 = sum(d4)/n/(prod(d4.^(1/n)));
M4 = Q*diag(d4)*Q';
M4(abs(M4)< 1e-12) = 0;  % getrid of small roundoff errors
M4 = sparse((M4+M4')/2);


% Matrix 5: Normally distributed with even smaller std:
d5_0 = (lam1-lamn)/100*randn(n-2,1)+(lam1+lamn)/2;
d5_0 = d5_0 + (meaneigs1 -mean(d5_0));

%d2_0 = d2_0 + m1-(lam1+lamn)/n;
d5 = [lamn;d5_0;lam1];
if min(d5)<lamn || max(d5)>lam1
    disp('Error with eigenvalues!!')
    keyboard
end
m5 = mean(d5);
std5 = std(d5);
omega5 = sum(d5)/n/(prod(d5.^(1/n)));
M5 = Q*diag(d5)*Q';
M5(abs(M5)< 1e-12) = 0;  % getrid of small roundoff errors
M5 = sparse((M5+M5')/2);


tol = 1e-6;
maxit = 2000;

xsol = randn(n,1);
temp = randn(n,1);
b = M1*xsol+M2*xsol+M3*xsol+M4*xsol+M5*xsol + 1e-7*temp/norm(temp);


[x1,flag1,relres1,iter1] = pcg(M1,b,tol,maxit);
[x2,flag2,relres2,iter2] = pcg(M2,b,tol,maxit);
[x3,flag3,relres3,iter3] = pcg(M3,b,tol,maxit);
[x4,flag4,relres4,iter4] = pcg(M4,b,tol,maxit);
[x5,flag5,relres5,iter5] = pcg(M5,b,tol,maxit);




fprintf('kappa always %g \n',kappa)
fprintf('Uniformly distributed eigs. ')
fprintf('Omega: %g Mean: %g  Std: %g pcg iter: %i relres: %e \n', omega1, m1,std1, iter1, relres1)

fprintf('Normally distributed eigs. ')
fprintf('Omega: %g Mean: %g Std: %g pcg iter: %i relres: %e \n', omega2, m2, std2, iter2, relres2)


fprintf('Normally distributed eigs. ')
fprintf('Omega: %g Mean: %g Std: %g pcg iter: %i relres: %e \n', omega3, m3, std3, iter3, relres3)

fprintf('Normally distributed eigs. ')
fprintf('Omega: %g Mean: %g Std: %g pcg iter: %i relres: %e \n', omega4, m4, std4, iter4, relres4)

fprintf('Normally distributed eigs. ')
fprintf('Omega: %g Mean: %g Std: %g pcg iter: %i relres: %e \n', omega5, m5, std5, iter5, relres5)

% [nRows,nCols] = size(base);
% [~,idx] = sort(rand(nRows,nCols),2);
% 
% %# convert column indices into linear indices
% idx = (idx-1)*nRows + ndgrid(1:nRows,1:nCols);
% 
% %# rearrange A
% B = base;
% B(:) = B(idx);

perc_o = @(x) (100 - ((x-1)*100)/(omega1-1));
perc_it = @(x) (100 - (x*100/iter1));

% For generating the table:

namefile = 'table_distributionofeigs.tex';

fmt0 = ' %s & %1.3e &  %1.4f & %i & - & - & %1.3e';
fmt = ' %s & %1.3e &  %1.4f & %i & %2.2f & 2.2f & %1.3e';

fid = fopen(namefile, 'w');

fprintf(fid,'%s\n', '\begin{tabular}{|l|r|r|r|r|r|r|} \hline');

fprintf(fid, '%s\n', ['Distribution & St. Deviation & $\omega$ & \# Iterations &' ...
    '\% Red. $\omega$ & \% Red. It. & Rel. Residual \\']);


fprintf(fid, ...
    '\n  \\hline \n');

fprintf(fid, fmt0, 'Uniform', std1, omega1, iter1, relres1);

fprintf(fid,'%s', '\\');

fprintf(fid, '\n');

fprintf(fid, fmt, 'Normal', std2, omega2, iter2, perc_o(omega2), perc_it(iter2), relres2);

fprintf(fid,'%s', '\\');

fprintf(fid, '\n');


fprintf(fid, fmt, 'Normal', std3, omega3, iter3, perc_o(omega3), perc_it(iter3), relres3);

fprintf(fid,'%s', '\\');

fprintf(fid, '\n');

fprintf(fid, fmt, 'Normal', std4, omega4, iter4, perc_o(omega4), perc_it(iter4), relres4);

fprintf(fid,'%s', '\\');


fprintf(fid, '\n');

fprintf(fid, fmt, 'Normal', std5, omega5, iter5, perc_o(omega5), perc_it(iter5), relres5);

fprintf(fid,'%s', '\\');


fprintf(fid, '\n');

fprintf(fid,' \\hline \n ');


fprintf(fid, '\\end{tabular}\n');

