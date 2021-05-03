function z_new=ACHRStep(prob0,K,z_old,old_points,center,is_warmup,solver)
%This function performs one step of ACHR sampling
%prob0 is the original LP problem
%z=K*x is sampled, generating z_new from z_old
%old_points is a matrix storing all previously sampled solutions
%center is the vector storing the center of all previously sampled
%solutions
%is_warmup == true: warm-up phase; false: main phase
%solver is the name of the algorithm used in solving the LP problem
if is_warmup == true
    deltav=rand_unit_vec(size(K,1));
else
    deltav=old_points(randi(size(old_points,1)),:)'-center(:);
    deltav=deltav/norm(deltav);
end
if nargin == 7
    drange=getACHRStepRange(prob0,K,z_old,deltav,solver);
else
    drange=getCHRRStepRange(prob0,K,z_old,deltav)
end
delta=rand()*(drange(2)-drange(1))+drange(1);
z_new=z_old+delta*deltav;
end