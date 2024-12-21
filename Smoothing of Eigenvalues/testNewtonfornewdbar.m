function d = testNewtonfornewdbar(A) 
%%% solve for new dbar  F(dbar) = Diag(dbar)B dbar - e = 0
%%   Here B =  Ainv.*Ainv;
%% INPUT A  positive definite
%% output  is d = diag sqrt Dbar 
%%     where Dbar = Diag(dbar)  
%%          min  omega (D Ainv D D Ainv D)
%%          Dinv * A * Dinv   is best scaline for Ax=b, A pos def
%%clear
%%rngsave = rng('shuffle');
%%n = 60;
%%
%%%   set up data
%%[q,~]=qr(randn(n));
%%Aeig = 10*(rand(n,1)+.1);
%%A = q*(diag(Aeig))*q';
%%A = (A + A')/2;
%%Ainv = inv(A);
%%B = Ainv.*Ainv;

%% initialize
n = length(A);
Ainv = inv(A);
B = Ainv.*Ainv;

F = @(d)( (diag(d)*B*d)-1 );
Jd = @(delta)( (diag(d)*B*delta)+(diag(delta)*B*d) );
In = eye(n);
omega = @(A)( (trace(A)/n)/(det_rootn(A)) );
%d = ones(n,1);
d = 1./sqrt(diag(B));
alpha = d'*B*d;
d = sqrt(n/alpha)*d;  % heuristic start pt
Fd = F(d);
resid = norm(Fd);
tol = 1e-10;
normB = norm(B,'fro');

iter = 0;
%while resid > tol*normB && iter <= 100
while resid > tol && iter <= 100
	iter = iter + 1;
	%J = zeros(n);
	%for ii = 1:n
	%	J(:,ii) = diag(d)*B*In(:,ii) + diag(In(:,ii))*B*d;
	%end
	J = diag(d)*B + diag(B*d);
	%if norm(J-Jrep,'fro') > 1e-8
	%	fprintf('temp test???  error J rep \n')
	%	norm(J-Jrep,'fro')
	%	keyboard
	%end
	dir = -J\Fd;
	d = d + dir;
	Fd = F(d);
	resid = norm(Fd);
end
%dapprox = sqrt(1./diag(B));
fprintf('iter = %i;  norm(F(d)) = %g;  d''Bd = %g; n = %i\n', ... 
	 iter,resid,d'*B*d,n) 
d = sqrt(d);  % for optimal output
