%%% Function for computing the w-optimal diagonal + t preconditioner
%%% in the form o f D_{+t} in Section 2.1.3 of the paper
%
% Input:
%       - W pos. def. matrix
%       - k - size of the triangular block
%
% Output:
%       - D = optimal preconditioner minimizing omega(D'*W*D)

function D = block_trir_preconditioner(W,k)
    
    n = length(W);
    
    tempW = W(1:k,1:k);
    R = chol(tempW);

    tempW = W(k+1:n,k+1:n);
    tempD = diag(diag(tempW).^(-1/2));
    
    D = blkdiag(inv(R),tempD);
    % D = zeros(n);
    % D(1:k,1:k) = inv(R);
    % 
    % 
    % D(k+1:n,k+1:n) = tempD;
end

