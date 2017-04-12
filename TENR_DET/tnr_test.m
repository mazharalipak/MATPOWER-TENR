
function [Tole,l,Iter]=tnr_test(datain,step)
%% Input Test case Structure and Intitial start.........................

testcase = datain;

Ldfactor1=1;
lambdan=0.0;
testcase.bus(:,3)=testcase.bus(:,3)*(1+lambdan*Ldfactor1);
testcase.bus(:,4)=testcase.bus(:,4)*(1+lambdan*Ldfactor1);
testcase.gen(:,2)=testcase.gen(:,2)*(1+lambdan*Ldfactor1);



mpc = loadcase(testcase);

mpc = ext2int(mpc);
sol = runpf(mpc);
sol=ext2int(sol);

%[ tnr ] = tnr_init(sol);

%% Input Data..........................

tnr = tnr_init(sol);
x = V2x_pf(tnr.V0, tnr);
l = 1;
z = [];

%% Newtpn Iterations ................................
tic;

Tole = 1;  
Iter = 1;
counter = 0;
while (Tole > 1e-2)
[ dx, dl,dMM,J] = tnr_step( x, l, z, tnr);
x=x+step*dx;
l=l+step*dl;

[~,a,~]=svds(J,1,0);

 Iter = Iter + 1;
 
      Tole=max(abs([dMM(1:end-1);a;dx;dl]));  
    counter = counter + 1;
    if counter ==250;
        break;
    end    
end
toc;


end


