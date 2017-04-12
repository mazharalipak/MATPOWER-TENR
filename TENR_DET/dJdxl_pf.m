function [ dJ,Gaa, Gav, Gva, Gvv] = dJdxl_pf( x, l, u, v, tnr, J)
%pf_dJ evaluate u'*dJ*v: sensitivity of the Jacobian
%   Detailed explanation goes here
    
    [pv, pq, npv, npq]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);
    
%     [J]=J_pf( x, l, tnr );   %% active this only for Determinant
%     Case......
    
    V = x2V_pf(x, tnr);
    
    lam = conj(F2S_pf(u, tnr));
    
    [Gaa, Gav, Gva, Gvv] = d2Sbus_dV2(tnr.Ybus, V, lam);

    dVa = v(1:npv+npq);
    dVm_pq = v(npv+npq+1:npv+2*npq);
    
    
    %% Determinant Case........
    
    dJJ= real ([ Gaa([pv;pq],[pv;pq])  Gav([pv;pq],pq); Gva(pq,[pv;pq])  Gvv(pq,pq)]);
    
    adjugate=adjugat(full(J));
    
    dJ=diag(sparse(adjugate)*dJJ);
    
    dJ=[dJ' 0];
    
    
    
%     dJ = real([
%         Gaa([pv;pq],[pv;pq])*dVa + Gav([pv;pq],pq)*dVm_pq; 
%         Gva(pq,[pv;pq])*dVa + Gvv(pq,pq)*dVm_pq;
%         0
%     ]');

end

