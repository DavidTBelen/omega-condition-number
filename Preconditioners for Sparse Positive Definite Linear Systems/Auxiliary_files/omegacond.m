%%% Function to compute the omega condition number based on the cholesky
%   decomposition:
%
% Input:
%       - A: positive definite matrix
%
% Output:
%       - wcond: omega condition number

function wcond = omegacond(A)
    
    n = length(A);
    
    f = trace(A)/n;
    R = chol(A);

    wcond = f/prod(diag(R).^(2/n));

end
