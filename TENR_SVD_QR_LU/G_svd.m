function [ G, u, v, E_M, F_M ] = G_svd(J, z)
%tr_svd Calculate the value of G and its derivatives
%   x - variable 
%   J - Jacobian function: evaluates J(x)

 [u, G, v] = svds(J, 1,0);
 
 E_M=sparse(eye(size(J)));
 F_M=E_M;
end

