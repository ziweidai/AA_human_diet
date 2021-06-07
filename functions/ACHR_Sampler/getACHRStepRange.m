function Range=getACHRStepRange(prob0,K,z,deltav,solver)
%Calculate range of 'jump' (perturbation) in one step of hit-and-run sampling
%prob0 stores the information of the original LP problem
%The transformed variable z=K*x is sampled, in which x subjects to
%constraints in prob0. prob0 = {a,blx,bux,blc,buc}
% blx <= x <= bux,
% blc <= a*x <= buc
% ivar is the coordinate to be perturbed
% z is the current value of z=K*x
% deltav is the direction for perturbation

A_Add=[K -deltav];
%epsilon=0;
prob1.a=[prob0.a zeros(size(prob0.a,1),1);A_Add];
prob1.blx=[prob0.blx;-1e8];
prob1.bux=[prob0.bux;1e8];
prob1.blc=[prob0.blc;z];
prob1.buc=[prob0.buc;z];
prob1.c=zeros(size(prob1.blx));
prob1.c(end)=1;

if nargin==5
    param.MSK_IPAR_OPTIMIZER=solver;
else
    param=[]
end
[~,res_min]=mosekopt('minimize echo(0)',prob1,param);
[~,res_max]=mosekopt('maximize echo(0)',prob1,param);
xmin=res_min.sol.bas.xx;
xmax=res_max.sol.bas.xx;
Range=[xmin(end) xmax(end)];
end