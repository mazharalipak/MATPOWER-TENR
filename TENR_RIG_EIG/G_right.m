function [ G, u, v ] = G_right(J, z)
%G_right Transversality condition based on right eigenvalue
%   Detailed explanation goes here
    n = size(J,1);

    G = [J*z; z'*z - 1];
    v = z;
    
    if issparse(J)
        u = sparse(1:n,1:n,ones(n,1), n, n+1);
    else
        u = [eye(n) zeros(n,1)];
    end

end

