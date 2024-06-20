%%% Function for computing the $\omega$-optimal incomplete upper triangular
%%% preconditioner
%
% Input:
%       - W pos. def. matrix
%       - k - size of the triangular block
%
% Output:
%       - D = optimal preconditioner minimizing omega(D'*W*D)

function D = i_upper_tri_preconditioner(W,k)
    
    n = length(W);
    
    tempW = W(1:k,1:k);
    R = chol(tempW);

    tempW = W(k+1:n,k+1:n);
    tempD = diag(diag(tempW).^(-1/2));
    
    D = blkdiag(inv(R),tempD);
    
end

