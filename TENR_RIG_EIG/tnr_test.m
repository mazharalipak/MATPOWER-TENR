
function [Tole,l,Iter]=tnr_test(datain,step)
%% Input Test case Structure and Intitial start.........................
testcase = datain;
mpc = loadcase(testcase);
mpc = ext2int(mpc);
sol = runpf(mpc);
sol=ext2int(sol);

%% Input Data..........................

tnr = tnr_init(sol);
[~, ~, npv, npq]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);

x = V2x_pf(tnr.V0, tnr);
l = 1;
z = 0.5*(ones((npv+2*npq),1));

%% Newtpn Iterations ................................
Tole = 1;  
Iter = 1;
counter = 0;
tic;

while (Tole > 1e-1)  
[ dx, dl,dz,dMM] = tnr_step( x, l, z, tnr);
x=x+step*dx;
l=l+step*dl;
z=z+step*dz;


 Iter = Iter + 1;
      Tole=max(abs([dMM;dx;dl]));  
    counter = counter + 1;
    if counter ==1;
        break;
    end    
end
toc;

end


