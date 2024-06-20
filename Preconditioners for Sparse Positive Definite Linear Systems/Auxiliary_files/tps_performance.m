%%% Given a set of variables 


function [tps,n_fails] = tps_performance(measure,fails)

nprec = size(measure,1);

%A = zeros(nprec,length(measure(1)));
for ii=1:nprec
    tps(:,ii) = cell2mat(measure(ii));
    tps_fails(:,ii) = cell2mat(fails(ii));
end

tps(tps_fails==1) = NaN;

n_fails = sum(tps_fails,1);

end