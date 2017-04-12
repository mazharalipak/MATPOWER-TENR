function [ G, u, v ] = G_det(J, z)
%tr_svd Calculate the value of G and its derivatives
%   x - variable 
%   J - Jacobian function: evaluates J(x)

 v= ones (size(J,1),1);
 u= ones (size(J,1),1);
 G=det(J);
end

