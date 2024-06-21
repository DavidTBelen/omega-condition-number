%%%% This generates a table showing number of condition number, #
%%%% iterations, and precision obtained by performig an iterative methd
%%%%for each one of the options for gamma

%namefile = ['../../../condnumbandQPlatexfiles.d/tableomega', solver, '.tex'];
%namefile = ['../../../condnumbandQPlatexfiles.d/', tablename, solver, '.tex'];
namefile = [tablename, solver, '.tex'];

fid = fopen(namefile, 'w');

fprintf(fid,'%s\n', '\begin{tabular}{|c|c|c|c|c|c|c|c|c|} \hline');

fprintf(fid, '%s', ...
     '$n$ & $\gamma$ & $\kappa(A(\gamma))$ & $\omega(A(\gamma))$ & Rel. Residual & No. Iterations & Time & T. $\gamma^{\ast}$ Spec & T. $\gamma^{\ast}$ Chol \\');

fprintf(fid, ...
    '\n  \\hline \n');


fmt1 = '%1.4e & %1.4e & %1.4e & %4.2f &  %2.4f & %2.4f & %2.4f';
fmt0 = '%1.4e & %1.4e & %1.4e & %4.2f &  %2.4f & - & -';

for nn = 1:length(sizesn)
n = sizesn(nn);
%row for gamma=0
fprintf(fid, '\\multirow{4}{*}{%d} & 0 &', n);

fprintf(fid, fmt0,...
    resultsn(nn,1,1), resultsn(nn,1,5),...
    resultsn(nn,1,3),resultsn(nn,1,2),resultsn(nn,1,4));

fprintf(fid, '%s \n \\cline{2-9} \n','\\');

%row for gamma=1
fprintf(fid, '& e &');

fprintf(fid, fmt0,...
    resultsn(nn,2,1), resultsn(nn,2,5),...
    resultsn(nn,2,3), resultsn(nn,2,2), resultsn(nn,2,4));

fprintf(fid, '%s \n \\cline{2-9} \n','\\');


%row for gamma=n
fprintf(fid, '& $u^{-2}$ &');

fprintf(fid, fmt0,...
    resultsn(nn,3,1), resultsn(nn,3,5),...
    resultsn(nn,3,3), resultsn(nn,3,2), resultsn(nn,3,4));

fprintf(fid, '%s \n \\cline{2-9} \n','\\');

%row for gammastar
fprintf(fid, '& $\\gamma^{\\ast}_p$ &');


fprintf(fid,fmt1,...
    resultsn(nn,4,1), resultsn(nn,4,5),...
    resultsn(nn,4,3), resultsn(nn,4,2), resultsn(nn,4,4),...
    timegammastarn(nn,1), timegammastarn(nn,2));

fprintf(fid,'%s \n \\hline \n ','\\');

end

fprintf(fid, '\\end{tabular}\n');
