function [ G, u, v, E_M] = G_lu(J, z) 
%tr_lu Calculate the value of G and its derivatives
%   x - variable 
%   J - Jacobian function: evaluates J(x)
% one can use option 1 or 2. The first one is an implementation based on a personally written LU function and might be a little slow.
% 2- second option is based on matlab lu decomposition and works just fine,....
% The code can be optimized if an impartial LU can be used for computing only smallest Pibot in U matrix

tol=1e-10;
[L_Mat, U_Mat, P_Mat, Q_Mat]=lucp(J,tol,'sparse');      %% option 1

%%

% thresh=1e-10;
% [L_Mat, U_Mat, P_Mat, Q_Mat]=lu(J, thresh, 'vector');  %% option 2

%%

New_UMat=bsxfun(@rdivide, U_Mat(1:end,:), diag(U_Mat));

[aa,bb]=min(abs(diag(U_Mat)));


b=zeros(size(J,1),1);                %% Right hand side ....(nth column of Identity matrix)
b(bb,1)=1;

v=(New_UMat)\b;                    

u=((((L_Mat))'))\b;

G=aa;

E_M=sparse(eye(size(J)));

end

