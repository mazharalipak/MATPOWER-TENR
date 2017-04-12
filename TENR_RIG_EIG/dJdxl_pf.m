function [ dJ,Gaa, Gav, Gva, Gvv] = dJdxl_pf( x, l, u, v, tnr)
%pf_dJ evaluate u'*dJ*v: sensitivity of the Jacobian
%   Detailed explanation goes here
    
    [pv, pq, npv, npq]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);
    
    V = x2V_pf(x, tnr);
    
    [a,b]=size(u);
    
    lam=zeros(npq+npv+1,b);
    
    for i=1:b;
    lam(:,i) = lam(:,i) + conj(F2S_pf(u(:,i), tnr));
    end
    
    dVa = v(1:npv+npq);
    dVm_pq = v(npv+npq+1:npv+2*npq);
        
    for i=1:b
    
    [Gaa, Gav, Gva, Gvv] = d2Sbus_dV2(tnr.Ybus, V, lam(:,i));
    dJ(:,i) = real([
        Gaa([pv;pq],[pv;pq])*dVa + Gav([pv;pq],pq)*dVm_pq; 
        Gva(pq,[pv;pq])*dVa + Gvv(pq,pq)*dVm_pq;
        0
    ]');

    end

end

