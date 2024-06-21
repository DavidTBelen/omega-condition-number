%%% This file generates a profile plot for the output of a chosen
%%% experiment

clear 
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));


% Chose table type: time, iterations, residuals, time_prec:

type = 'time_prec';


% Chose the data set:
which_experiment = 'selected_matrices';

% Save the plot:
save_table = true;

load_data = true;

if load_data == true
    dataname = ['output_exp_', which_experiment];
    load(dataname)
    
    data = load([which_experiment,]); %array with name of matrices
    matrices_names = data.selected_matrices_list;


    
    data = output;
else

end

% Load fails:

n = output.nproblems;
sizes = output.problem_sizes;
nnzW = output.problem_nnz;

fails = zeros(n,4);
fails(:,1) = output.no_prec.fails;
fails(:,2) = output.diag_prec.fails;
fails(:,3) = output.dbt_prec.fails;
fails(:,4) = output.ichol_prec.fails;

if strcmp(type, 'time')
    measure = zeros(n,4);
    measure(:,1) = output.no_prec.time;
    measure(:,2) = output.diag_prec.time;
    measure(:,3) = output.dbt_prec.time;
    measure(:,4) = output.ichol_prec.time;

    % namefile = '..\..\Latexrevisingasofnov9COAPSubmission/table_time.tex';
    namefile = 'table_time.tex';

    
    fmt0 = '%i & %i ';  %1.4e & %4.2f &  %2.4f & %2.4f & %2.4f';
    fmt1 = ' & %.2f ';
    fmt2 = ' & $>$%.2f ';
   
elseif strcmp(type, 'iterations')
    measure = zeros(n,4);
    measure(:,1) = output.no_prec.it;
    measure(:,2) = output.diag_prec.it;
    measure(:,3) = output.dbt_prec.it;
    measure(:,4) = output.ichol_prec.it;

    % namefile = '..\..\Latexrevisingasofnov9COAPSubmission/table_iterations.tex';
    namefile = 'table_iterations.tex';

    
    fmt0 = '%i & %i ';  %1.4e & %4.2f &  %2.4f & %2.4f & %2.4f';
    fmt1 = ' & %.1i ';
    fmt2 = ' & $>$%.1i ';
    fmt3 = '& -';

elseif strcmp(type, 'residuals')
    measure = zeros(n,4);
    measure(:,1) = output.no_prec.res;
    measure(:,2) = output.diag_prec.res;
    measure(:,3) = output.dbt_prec.res;
    measure(:,4) = output.ichol_prec.res;

    % namefile = '..\..\Latexrevisingasofnov9COAPSubmission/table_residuals.tex';
    namefile = 'table_residuals.tex';
    
    fmt0 = '%i & %i ';  %1.4e & %4.2f &  %2.4f & %2.4f & %2.4f';
    fmt1 = ' & %1.3e ';
    fmt2 = ' & - ';

elseif strcmp(type, 'time_prec')
    measure = zeros(n,4);
    measure(:,1) = zeros(n,1);
    measure(:,2) = output.diag_prec.time_prec;
    measure(:,3) = output.dbt_prec.time_prec;
    measure(:,4) = output.ichol_prec.time_prec;

    % namefile = '..\..\Latexrevisingasofnov9COAPSubmission/table_time_prec.tex';
    namefile = 'table_time_prec.tex';

    
    fmt0 = '%i & %i ';  %1.4e & %4.2f &  %2.4f & %2.4f & %2.4f';
    fmt1 = ' & %1.3e ';
    fmt2 = ' & %1.3e ';
    fmt3 = '& %i';

end


fid = fopen(namefile, 'w');

if strcmp(type, 'time_prec')

    fprintf(fid,'%s\n', '\begin{tabular}{|l|rr|r|r|r|} \hline');
    
    fprintf(fid, '%s', ...
         'name & $n$ & $nnz(W)$  & DIAG & ITRIU & ICHOL \\');
else
    fprintf(fid,'%s\n', '\begin{tabular}{|l|rr|r|r|r|r|} \hline');
    
    fprintf(fid, '%s', ...
         'name & $n$ & $nnz(W)$ & NONE & DIAG & ITRIU & ICHOL \\');
end

fprintf(fid, ...
    '\n  \\hline \n');


for ii = 1:n

    fprintf(fid, '\\verb|%s| &', matrices_names(ii));

    fprintf(fid, fmt0,  sizes(ii),  nnzW(ii));
    

    for jj=1:4

        if strcmp(type, 'time_prec') &&  jj==1
        
        else

           
            if fails(ii,jj) == 1

                if strcmp(type, 'iterations') && measure(ii,jj) == 0
                    fprintf(fid,fmt3,  measure(ii,jj));
                else
                    fprintf(fid,fmt2,  measure(ii,jj));
                end
            else
                fprintf(fid,fmt1,  measure(ii,jj));
            end
        end
    end

    fprintf(fid, '%s \n', '\\');

end

fprintf(fid,' \\hline \n ');


fprintf(fid, '\\end{tabular}\n');




