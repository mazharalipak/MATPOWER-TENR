function [ G, u, v, E_M] = G_qr(J, z)
%tr_qr Calculate the value of G and its derivatives
%   x - variable 
%   J - Jacobian function: evaluates J(x)

[Q,R,E]=qr(J);

E_M=E;                           % Permutation Matix for 

b=zeros(size(J,1),1);             %% Right hand side ....(nth column of Identity matrix)
b(end,end)=1;

RR=R;
RR(end,end)=1;

v=sparse(RR\b);                   %% Finding .......(R_\tilde(:,END))^{-1}
 
u=Q(:,end);

G=R(end,end);


end

