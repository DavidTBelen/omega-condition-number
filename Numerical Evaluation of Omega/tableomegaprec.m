%%%% This generates a table with comparision of different methods for
%%%% obtaining omega condition number.

namefile = '../../condnumbandQPlatexfiles.d/tableomegaprec.tex';

fid = fopen(namefile, 'w');

fprintf(fid,'%s', '\begin{tabular}{|c|c|');

for ii=1:length(xkappa)
    fprintf(fid,'c|');
end

fprintf(fid, '} \n \\hline \n');


fprintf(fid, '%s', ...
     '$n$ & Fact.  ');

for ii=1:length(xkappa)
    fprintf(fid, '& order %s 1e%i', '$\kappa$', xkappa(ii));
end
fprintf(fid, '%s \n \\hline \n','\\');

fmt = [];
for ii=1:length(xkappa)
    fmt = [fmt, ' & %f '];
end


for jj = 1:length(n)
    nn = n(jj);
    %for eig
    fprintf(fid, '\\multirow{3}{*}{%d} & eig ', nn);
    
    for kk =1:length(xkappa)
        fprintf(fid, '& %1.4e', meanresultprec(jj,kk,1));
    end

    
    fprintf(fid, '%s \n \\cline{2-%i} \n','\\',length(xkappa)+2);
    
    
    
    % for chol 
    fprintf(fid, ' & R ');
    
       
    for kk =1:length(xkappa)
        fprintf(fid, '& %1.4e', meanresultprec(jj,kk,2));
    end

    
    fprintf(fid, '%s \n \\cline{2-%i} \n','\\',length(xkappa)+2);


    % for LU
    
    fprintf(fid, ' & LU ');
    
    for kk =1:length(xkappa)
        fprintf(fid, '& %1.4e', meanresultprec(jj,kk,3));
    

    end
    
    fprintf(fid,'%s \n \\hline \n ','\\');

end

fprintf(fid, '\\end{tabular}\n');
% % % 
% % % 
% % % 
% % % for nn = 1:length(sizesn)
% % % n = sizesn(nn);
% % % %row for gamma=0
% % % fprintf(fid, '\\multirow{3}{*}{%d} & 0 &', n);
% % % 
% % % fprintf(fid, fmt,...
% % %     resultsn(nn,1,1), resultsn(nn,1,5),...
% % %     resultsn(nn,1,3),resultsn(nn,1,2),resultsn(nn,1,4));
% % % 
% % % fprintf(fid, '%s \n \\cline{2-5} \n','\\');
% % % 
% % % %row for gamma=1
% % % fprintf(fid, '& e &');
% % % 
% % % fprintf(fid, fmt,...
% % %     resultsn(nn,2,1), resultsn(nn,2,5),...
% % %     resultsn(nn,2,3), resultsn(nn,2,2), resultsn(nn,2,4));
% % % 
% % % fprintf(fid, '%s \n \\cline{2-5} \n','\\');
% % % 
% % % %row for gammastar
% % % fprintf(fid, '& $\\gamma^{\\ast}$ &');
% % % 
% % % 
% % % fprintf(fid,fmt1,...
% % %     resultsn(nn,3,1), resultsn(nn,3,5),...
% % %     resultsn(nn,3,3), resultsn(nn,3,2), resultsn(nn,3,4),...
% % %     timegammastarn(nn,1), timegammastarn(nn,2));
% % % 
% % % fprintf(fid,'%s \n \\hline \n ','\\');
% % % 
% % % end
% % % 
% % % fprintf(fid, '\\end{tabular}\n');
