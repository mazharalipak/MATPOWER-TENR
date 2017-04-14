
function [Tole,l,Iter]=tnr_test(datain,step,Lambda)
%% Input Test case Structure and Intitial start.........................
t=1;                              % Slack variable for \lambda >0 .....................

%% Initial seed / germ for starting TENR Algorithm...........................
testcase = datain;
testcase.bus(:,3)=testcase.bus(:,3)*(1+Lambda);
testcase.bus(:,4)=testcase.bus(:,4)*(1+Lambda);
testcase.gen(:,2)=testcase.gen(:,2)*(1+Lambda);

%% Loading Data structure for computing solutions..............

mpc = loadcase(testcase);
mpc = ext2int(mpc);
sol = runpf(mpc);
sol=ext2int(sol);

%% Input Data..........................

tnr = tnr_init(sol);
type=sol.bus(:,2);
sol.bus(:,2)=type;
x = V2x_pf(tnr.V0, tnr);
l = 1;
z = [];


%% Reactive generation bound on generator buses.....

vec_1=sol.gen(:,1);

Qg_max=zeros(tnr.npq+tnr.npv+1, 1);
Qg_max(vec_1,1)=(sol.gen(:,4)./100);

Qg_min=zeros(tnr.npq+tnr.npv+1,1);
Qg_min(vec_1,1)=(sol.gen(:,5)./100);

%% Newtpn Iterations ................................
tic;

Tole = 1;  
Iter = 1;
counter = 0;
while (Tole > 1e-2)  
 
%% TENR Iteration 

[ V ]= x2V_pf( x, tnr);

%V_abs=abs(V);

Qgen_cal= imag(V .*conj(tnr.Ybus * V))+ (sol.bus(:,4)./100);

%%

   if Iter <= 150 && Iter > 1    % Only checked up to 7th iterations..
        for n = tnr.pv
                 %Qgen_cal=Qgen_cal(n);
                if Qgen_cal(n) < Qg_min(n)
                    %tnr.V0(n) = tnr.V0(n) + 0.211;
                    type(n)=1;
                elseif Qgen_cal(n) > Qg_max(n)
                    type(n)=1;
                    %tnr.V0(n) = tnr.V0(n) - 0.211;
                end
        end
   end
   
%%
[ dx, dl,dt,dMM,~] = tnr_step( x, l,z,t, tnr,sol);

%% Updating state variables......

x=x+step*dx;
l=l+step*dl;
t=t+step*dt;
 
%% Checking Tolerance.......

 Iter = Iter + 1;
      Tole=max(abs([dMM;dx;dl]));  
    counter = counter + 1;
    if counter ==100;
        break;
    end    
end

toc;
end



