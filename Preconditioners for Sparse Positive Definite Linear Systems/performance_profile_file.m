%%% This file generates a profile plot for the output of a chosen
%%% experiment

clear all
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));


% Chose the data set:
which_experiment = 'selected_matrices';

% Save the plot:
save_plot = true;

load_data = true;

if load_data == true
    dataname = ['output_exp_', which_experiment];
    load(dataname)
    
    data = output;
else

end



%%%%%%%%%%%%%%%%%%%%%

% Define the data for the performance profile:

fails = cell(4,1);
fails(1) = {output.no_prec.fails};
fails(2) = {output.diag_prec.fails};
fails(3) = {output.dbt_prec.fails};
fails(4) = {output.ichol_prec.fails};


measure_time = cell(4,1);
measure_time(1) = {output.no_prec.time};
measure_time(2) = {output.diag_prec.time};
measure_time(3) = {output.dbt_prec.time};
measure_time(4) = {output.ichol_prec.time};


measure_it = cell(4,1);
measure_it(1) = {output.no_prec.it};
measure_it(2) = {output.diag_prec.it};
measure_it(3) = {output.dbt_prec.it};
measure_it(4) = {output.ichol_prec.it};



% Performance profile times

fig1 = figure(1);
clf
% figname  = ['..\..\Latexrevisingasofnov9COAPSubmission/perf_time_', which_experiment, '.png'];
figname  = ['perf_time_', which_experiment, '.png'];
[T, n_fails] = tps_performance(measure_time,fails);
perf(T,1,'Time',n_fails)


if save_plot == true
    exportgraphics(fig1,figname)
end

% Performance profile iterations:

fig2 = figure(2);
clf
% figname  = ['..\..\Latexrevisingasofnov9COAPSubmission/perf_it_', which_experiment, '.png'];
figname  = ['perf_it_', which_experiment, '.png'];
[T, n_fails] = tps_performance(measure_it,fails);
perf(T,1,'Iterations',n_fails)

if save_plot == true
    exportgraphics(fig2,figname)
end