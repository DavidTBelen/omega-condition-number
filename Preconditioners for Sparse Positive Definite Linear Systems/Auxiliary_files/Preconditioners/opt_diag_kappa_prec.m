function D = opt_diag_kappa_prec(W)
    
    n = length(W);
    %dhat = (1/n)*ones(n,1);  %starting vector
    dhat = rand(n,1);  %starting vector
    dhat = dhat/norm(dhat);
    dhat = norm(W,'fro')*dhat/2;
    S = sqrtm(full(W));
    Wd = @(d)( (S*diag(d)*S) );   % Wdiag(d) = Sdiag(d)S for tr/det/eig fns

    %gomega = @(d)(  det_rootn(diag(d))  );   % det W constant
    fkappa = @(d)( ( lambda_max(Wd(d))  )  );
    gkappa =  @(d)( (lambda_min(Wd(d)))  );
    Fkappa = @(d)(  fkappa(d)/gkappa(d) );
    %qhatkappa = Fkappa(dhat);  % initial data
    
    f = fkappa;
    g = gkappa;
    CondF = Fkappa;
    %qhat = qhatkappa; %not used???????????????
    
    
    %%%%%%%% find the optimal q using Dinkelbach
    dsnew = dhat;
    Fq = -inf;    % minimization so F is increasing and is negative
    toler = 1e-12;
    maxiter = 10;
    iter = 1;
    Fqs = [];
    qs = [];
    reltoler =  toler*(norm(W,'fro'));
    qs(iter) = CondF(dsnew) + 10;   % ensure f(dsnew) - qs*gsnew < 0 as f>0,g>0
    M = 10;  % to start negative
    Fqs(iter) =  f(dsnew) - qs(iter)*g(dsnew) - M;
    
    tic
    while Fq  < -reltoler && iter <= maxiter
    
	    cvx_begin
	    cvx_precision best
	    variable ds(n,1) nonnegative
	    dual variable lam
	    minimize( f(ds) - qs(iter)*g(ds));
	    subject to
		    lam : norm(ds) <= norm(W,'fro');
	    cvx_end
    
    
	    dsnew = ds;
	    Fq = cvx_optval;
	    %figure(2)
	    %clf
	    %% use -Fqs for semilogy that needs Pos els
	    %plot(qs(1:iter),Fqs(1:iter),'-x')
	    %semilogy(qs(1:iter),-min(0,Fqs(1:iter)),'-x')
	    %legend('q value','Fq value','location','best')
	    %xlabel('q values')
	    %ylabel('-F(q) values')
	    %title(['q vs -Fq in ',num2str(iter), ...
		%    ' iters; Dinkelbach alg. for Xsp'])
	    %drawnow
	    iter = iter + 1;
	    Fqs(iter) = Fq;
	    qs(iter) = CondF(dsnew);
    end
    D = ds;
end