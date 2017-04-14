function [ dx, dl,dt,dMM,J] = tnr_step( x, l,z,t, tnr,sol)
%tnr_step Summary of this function goes here
%   Detailed explanation goes here


%%  Routines for Creating Elements of Extended Jacobain.....

[nx, ~] = deal(length(x), length(z));
    
F = tnr.F(x, l, tnr);
[J,JJ] = tnr.J(x, l, tnr,sol);

[ G, u, v, E_M] = tnr.G(J, z);
dGdxl = tnr.dJdxl(x, l, u, v, E_M, tnr);
        
%% enforcing Lambda>0 condition.............

Delta_M=l-t^2;                                     % nonlinear equality equation for Lambda >0 ................
        
ra=zeros(1,size(JJ,2)+1);
ra(1,end-1)=1;
ra(1,end)=-2*t;
        
%% Extended Jacobian Matrix................
        
fullJJ=[JJ zeros(size(JJ,1),1); dGdxl 0; ra];
                     
%% Solving Linear system of Equation Newton Step....

dMM=[F;G;Delta_M];
        
dv = - fullJJ \ dMM;
    
dx = dv(1:nx);
dl = dv(nx+1);
dt = dv(end,1);
end

