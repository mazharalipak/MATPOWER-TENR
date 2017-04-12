function [ S ] = F2S_pf( f, tnr )
%pf_F2S Converts the F vector to powers (ignoring the reactive power on
%generators)
%   Detailed explanation goes here
    [pv, pq, npv, npq]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);

    S = zeros(size(tnr.Ybus,1),1);
    S(pv) = f(1:npv);
    S(pq) = f(npv+1:npv+npq) + 1i*f(npv+npq+1:npv+2*npq);

end

