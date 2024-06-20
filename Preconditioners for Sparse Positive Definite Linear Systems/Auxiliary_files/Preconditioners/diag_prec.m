%%% Function for computing the w-optimal diagonal preconditioner
%
% Input:
%       - W pos. def. matrix
%
% Output:
%       - D = Diag preconditioner minimizing omega(D'*W*D)

function D = diag_prec(W)
    d = diag(W);
    D = diag(d.^(-1/2));
end

