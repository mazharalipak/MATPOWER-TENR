function [J,JJ] = J_pf( x, l, tnr )
%J_pf Calculate the residual of power flow Jacobian
%   Detailed explanation goes here

    [pv, pq]  = deal(tnr.pv, tnr.pq);

    V = x2V_pf(x, tnr);
    
    [dSbus_dVm, dSbus_dVa] = dSbus_dV(tnr.Ybus, V);
    
    j11 = real(dSbus_dVa([pv; pq], [pv; pq]));
    j12 = real(dSbus_dVm([pv; pq], pq));
    j21 = imag(dSbus_dVa(pq, [pv; pq]));
    j22 = imag(dSbus_dVm(pq, pq));
        
        J=[  j11 j12 ;
            j21 j22 ];
        
        vec_F=[-real(tnr.dS([pv;pq]));-imag(tnr.dS(pq))];
        
        JJ= [ J zeros(size(J)) vec_F];

%         JJ=[j11 j12 zeros(size(j11)) -real(tnr.dS([pv;pq]));
%             j21 j22 zeros(size(j21)) -imag(tnr.dS(pq))];
        
end

