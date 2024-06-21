function tps = performanceprofilegen(results,timegammastar,optionsperf)
%%% Given a vector of results of the lsqr experiment, generates a dot plot
%    results = results(length(sizen),np,i,j)
%    i = 1 -> gamma = 0; 2 -> gamma = 1; 3 -> gamma^star
%    j = 1 -> condest; 2-> iter; 3 -> relres; 4 -> time; 5 -> omega cond
%%%%%%%%%%%%%%%%%%%%%
% options.lsqr or options.cgs
% optionsperf.maxit -- max iter employed

ndim = size(results,1);
nprob = size(results,2);
%maxit = optionsperf.maxit;
tol = optionsperf.tol;



tps = zeros(ndim*nprob,4);


for ii=1:ndim
    for jj=1:nprob
        for kk=1:3
            if results(ii,jj,kk,3) < tol
                tps((ii-1)*nprob +jj,kk) =  results(ii,jj,kk,4);
            else
                tps((ii-1)*nprob+jj,kk) = NaN;
            end
        end
        if results(ii,jj,4,3) < tol
            tps((ii-1)*nprob+jj,4) = results(ii,jj,4,4) + timegammastar(ii,jj,1);
        else
            tps((ii-1)*nprob+jj,4) = NaN;
        end
    end
end

end %function

%rps = NaN*ones(ndim*nprob,3); %3 because we just consider 3 choices of gamma

 
% tps = zeros(ndim,nprob,3);
% tps(:,:,1:2) = results(:,:,1:2,4);
% tps(:,:,3) = results(:,:,3,4) + timegammastar(:,:,2);
% 
% for ii=1:ndim
%     for jj=1:nprob
%     miniijj = min(tps(ii,jj,:));
%         for kk=1:3
%             if reesults(ii,jj,kk,3) > tol
%                 rps(ii*jj,kk) = tps(ii,jj,kk)/miniijj;
%             end
%         end
%     end
% end





