function perf(T,logplot,title4str,n_fails)
%PERF    Performace profiles
%
% PERF(T,logplot)-- produces a performace profile as described in
%   Benchmarking optimization software with performance profiles,
%   E.D. Dolan and J.J. More', 
%   Mathematical Programming, 91 (2002), 201--213.
% Each column of the matrix T defines the performance data for a solver.
% Failures on a given problem are represented by a NaN.
% The optional argument logplot is used to produce a 
% log (base 2) performance plot.
if (nargin < 2) logplot = 0; end

colors  = [ [1 0 0];[0 0 1];[0 1 0];[0 1 1];[1 0 1];
            [0 0 0.5];[0 0.5 0];[0 0.5 0.5];[0.5 0 0];[0.5 0 0.5];[0.5 0.5 0];
            [0.3 0.9 0.5];[0.7 0.3 0.5];[0.5 0.3 0.5];[0.3 0.3 0.7];[0.2 0.5 0.9];[0.1 0.4 0.3]
          ];
%lines   = cellstr(char( '-.', '--', ':', '-'));
LineStyles   = ["-.","--","-",":","-."];  %[ '-.' '--' ':' '-'];
%markers = ['o' 's'  'x' 'd' '>' 'v' '*'];

[np,ns] = size(T);

% Minimal performance per solver

minperf = min(T,[],2);

% Compute ratios and divide by smallest element in each row.

r = zeros(np,ns);
for p = 1: np
  r(p,:) = T(p,:)/minperf(p);
end

if (logplot) r = log2(r); end

max_ratio = max(max(r));

% Replace all NaN's with twice the max_ratio and sort.

topline = 2*max_ratio; % To add line on the top of plot
r(find(isnan(r))) = topline;
r = sort(r);


% Plot stair graphs with markers.
clf;
for s = 1: ns
 [xs,ys] = stairs([0; r(:,s); topline],[0; [1:np]'/np; 1]); 
 if sum(xs == topline) == length(xs)-1
     xs = linspace(0,topline,length(xs))';   
     ys = zeros(length(ys),1);
 end

     
 temp =  s;  %floor((s-1)./3)+1;
 %set(gca, 'LineStyle', '-');
 plot(xs,ys,  LineStyles(temp) ,'Color', colors(s,:),'LineWidth' ,1);
 hold on;
end

legend_labels_1 = "NONE (" +string(n_fails(1))+ " failures)";
legend_labels_2 = "DIAG (" +string(n_fails(2))+ " failures)";
legend_labels_3 = "ITRIU (" +string(n_fails(3))+ " failures)";
legend_labels_4 = "ICHOL(1) (" +string(n_fails(4))+ " failures)";
legend_labels_5 = "ICHOL(2) (" +string(n_fails(5))+ " failures)";


legend(legend_labels_1,	legend_labels_2, legend_labels_3, legend_labels_4,...
    legend_labels_5);
legend('location','southeast')
if (logplot) 
    xlabel(strcat('log_2(\tau)')); 
    ylabel('\rho_{\gamma}(\tau)');
else
    xlabel('\tau');
    ylabel('\rho_{\gamma}(\tau)')
end

fontSize = 14; 
caption = sprintf([title4str, ' performance profile']);
title(caption, 'FontSize', fontSize);

% Axis properties are set so that failures are not shown,
% but with the max_ratio data points shown. This highlights
% the "flatline" effect.


axis([ 0 1.1*max_ratio 0 1 ]);

% Legends and title should be added.