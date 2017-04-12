function [F] = F_pf( x, l, tnr )
%pfF Calculate the residual of power flow Jacobian
%   Detailed explanation goes here

    [pv, pq]  = deal(tnr.pv, tnr.pq);

    V = x2V_pf(x, tnr);
    
    residual =  V .* conj(tnr.Ybus * V) - tnr.S0 - l*tnr.dS;
    F = [ 
      real(residual(pv));
      real(residual(pq));
      imag(residual(pq));
    ];
end

