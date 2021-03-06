%-------------------------------------------------------------------------
% Sample AA compositions in different diets
%-------------------------------------------------------------------------

%% Define variables
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

%% Process the sampled values and show violin plots
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
    violinplot(aa_ratio_samp,DietNames,'ShowData',false);
    xtickangle(45);
    if i < 16
        xticklabels([]);
    end
    title(AANames{i});
    xlim([0 11]);
    box on;
end
mean_aa_ratio=zeros(10,18);
for i=1:10
    a=list_aa_samp{i};
    mean_aa_ratio(i,:)=mean(a./sum(a,2));
end

%% Show violin plots for AAs most variable across dietary patterns
[~,idx] = sort(F_ANOVA_AA_Diet,'descend');
figure;
for i=1:4
    aa_ratio_samp=zeros(nFinal,10);
    for j=1:10
        a=list_aa_samp{j};
        aa_ratio_samp(:,j)=a(:,idx(i))./sum(a,2);
    end
    subplot(2,2,i);
    violinplot(aa_ratio_samp,DietNames,'ShowData',false);
    xtickangle(45);
    title(AANames{idx(i)});
    xlim([0 11]);
    ylabel('Abundance [g/g total AA]');
    box on;
end

%% Plot average fractions of AAs in dietary patterns
data = (mean_aa_ratio-min(mean_aa_ratio))./(max(mean_aa_ratio)-min(mean_aa_ratio));
figure;

cmap_now = brewermap(100,'RdYlGn');
heatmap_cluster(data',AANames,DietNames,[0 1],cmap_now(75:-1:26,:));
colorbar('Ticks',[0 1],'TickLabels',{'Row min','Row max'});
title('Relative abundance of amino acids in diets');
clear data cmap_now

%% PCA of amino acid compositions of diets
diet_tag=reshape(repmat(DietNames,5000,1),50000,1);
aa_ratio_downsample=zeros(50000,18);
for i=1:10
    rp=randperm(50000);
    a=list_aa_samp{i};
    aa_ratio_downsample((i-1)*5000+1:i*5000,:)=a(rp(1:5000),:)./sum(a(rp(1:5000),:),2);    
end

[~,score,~,~,explained,~]=pca(aa_ratio_downsample);
figure;
color=brewermap(11,'Set3');
for i=1:10
    data=score((i-1)*5000+1:i*5000,:);
    scatter3(data(:,1),data(:,2),data(:,3),5,color(i,:),'filled');
    hold(gca,'on');
end
box on;
grid on;
xlabel(sprintf('PC1 (%.1f%%)',explained(1)));
ylabel(sprintf('PC2 (%.1f%%)',explained(2)));
zlabel(sprintf('PC3 (%.1f%%)',explained(3)));
legend(DietNames);
title('PCA of relative amino acid compositions of diets');

