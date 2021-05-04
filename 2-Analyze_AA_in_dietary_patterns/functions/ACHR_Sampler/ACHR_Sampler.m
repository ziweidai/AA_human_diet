function ZFinal=ACHR_Sampler(prob,K,z0,nWarmup,nFinal)

%ACHR sampling on the nutrient variables
%prob: data structure storing linear constraints
%K: coefficient matrix for the variables to be sampled: z=K*x is sampled
%instead of x
%z0: z0=K*x0 is the starting point of the sampling
%nWarmup: number of warm-up points
%nFinal: number of points to be sampled

solver='MSK_OPTIMIZER_FREE_SIMPLEX';

zdim=size(K,1);
ZSamp=zeros(nWarmup+nFinal,zdim);
z_old=z0(:);
ZSamp(1,:)=z_old';

for i=1:nWarmup-1
    z_new=ACHRStep(prob,K,z_old,ZSamp(1:i,:),zeros(1,zdim),true,solver);
    ZSamp(i+1,:)=z_new';
    z_old=z_new;
end

center=mean(ZSamp(1:nWarmup,:));
for i=1:nFinal
    z_new=ACHRStep(prob,K,z_old,ZSamp(1:nWarmup+i-1,:),center,false,solver);
    center=((nWarmup+i-1)*center+z_new')/(nWarmup+i);
    ZSamp(i+nWarmup,:)=z_new';
    z_old=z_new;
end

ZFinal=ZSamp(nWarmup+1:end,:);
end