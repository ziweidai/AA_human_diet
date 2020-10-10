%-------------------------------------------------------------------------
% Sample AA compositions in different diets
%-------------------------------------------------------------------------

%% Define variables
DefineDiets;
%Coefficients for sampling AA compositions
K_AA_all=FoodMatrix(CompletePos,AAPos)'/100;
%K_totalAA=sum(K_AA_all);
nWarmup=1000;
nFinal=50000;

%% Sample AA composition of diets
list_aa_samp=cell(1,10); %Used to be list_aa_ratio_samp which is obtained 
                         %from independently sample (AA,tAA) pairs
for i=1:10
    Diet=DietList{i};
    DietNames(i)
    prob_samp_diet.blx=zeros(2335,1);
    prob_samp_diet.bux=1000*ones(2335,1); %No more than 1kg for one food
    prob_samp_diet.blc=[Diet.ConstraintLBs;0];
    prob_samp_diet.buc=[Diet.ConstraintUBs;3000];  %No more than 3kg total food
    prob_samp_diet.a=[Diet.ConstraintMatrix(:,CompletePos);ones(1,2335)];
    if i<9
        prob_samp_diet.blc=[prob_samp_diet.blc;1800];
        prob_samp_diet.buc=[prob_samp_diet.buc;2200];
        prob_samp_diet.a=[prob_samp_diet.a;FoodMatrix(CompletePos,143)'/100];
    end
    x0=[];
    for j=1:100
        prob_samp_diet.c=rand(2335,1)-0.5;
        [~,res1]=mosekopt('minimize echo(0)',prob_samp_diet);
        [~,res2]=mosekopt('maximize echo(0)',prob_samp_diet);
        x0=[x0 (res1.sol.bas.xx+res2.sol.bas.xx)/2];       
    end
    x0=mean(x0')';
    clear res1 res2
    list_aa_samp{i}=ACHR_Sampler(prob_samp_diet,K_AA_all,K_AA_all*x0,nWarmup,nFinal);
end

%% Process the sampled values
F_ANOVA_AA_Diet=zeros(18,1);
figure;
for i=1:18
    aa_ratio_samp=zeros(nFinal,10);
    for j=1:10
        a=list_aa_samp{j};
        aa_ratio_samp(:,j)=a(:,i)./sum(a,2);
    end
    [~,tbl]=anova1(aa_ratio_samp,[],'off');
    F_ANOVA_AA_Diet(i)=tbl{2,5};
    subplot(6,3,i);
    [~,idx]=sort(median(aa_ratio_samp));
    violinplot(aa_ratio_samp(:,idx),DietNames(idx),'ShowData',false);
    xtickangle(45);
    title(AANames{i});
    xlim([0 11]);
    box on;
end
mean_aa_ratio=zeros(10,18);
for i=1:10
    a=list_aa_samp{i};
    mean_aa_ratio(i,:)=mean(a./sum(a,2));
end
figure;
plotFractions(mean_aa_ratio,DietNames,AANames);

%% Generate t-SNE plot of sampled AA compositions of diets
diet_tag=reshape(repmat(DietNames,5000,1),50000,1);
aa_ratio_downsample=zeros(50000,18);
for i=1:10
    rp=randperm(50000);
    a=list_aa_samp{i};
    aa_ratio_downsample((i-1)*5000+1:i*5000,:)=a(rp(1:5000),:)./sum(a(rp(1:5000),:),2);    
end
x=tsne(aa_ratio_downsample);
color=brewermap(11,'Set3');
figure;
gscatter(x(:,1),x(:,2),diet_tag,color,[],[]);

%% Plot AA composition of example diets under each dietary scheme
n_exp=3; %Number of exemplar diets in each diet type
for i=1:10
    x=aa_ratio_downsample((i-1)*5000+1:i*5000,:);
    d=squareform(pdist(x));
    exp_pos=find(max(d)==max(d(:)));
    for j=1:n_exp-2
        mind=min(d(exp_pos,:));
        [~,idx]=find(mind==max(mind));
        exp_pos=[exp_pos idx];
    end
    subplot(5,2,i);
    if i<=8
        plotFractions(x(exp_pos,:),{'Diet 1','Diet 2','Diet 3'},[]);
    else
        plotFractions(x(exp_pos,:),{'Diet 1','Diet 2','Diet 3'},AANames);
    end
    title(DietNames{i});
    yticklabels({'Diet 1','Diet 2','Diet 3'});
end