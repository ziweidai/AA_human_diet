%-------------------------------------------------------------------------
% Sample carbohydrate and lipid compositions in different diets
%-------------------------------------------------------------------------
DefineDiets;
SampVarNames={'Carbohydrate, by difference,g','Fiber, total dietary,g',...
    'Sugars, total,g','Fatty acids, total saturated,g',...
'Fatty acids, total polyunsaturated,g','Fatty acids, total monounsaturated,g'};
[~,pos]=ismember(SampVarNames,Nutrients);

%% Sample carb/fat composition in diets
%Coefficients for sampling carb and fat compositions
K_CarbFat_all=FoodMatrix(CompletePos,pos)'/100;
%K_totalAA=sum(K_AA_all);
nWarmup=1000;
nFinal=50000;
n_foods=length(CompletePos);

%----------------------------------------------------
%Calculate distance to other diets except for USDA
%----------------------------------------------------
list_carbfat_samp=cell(1,10); %Used to be list_aa_ratio_samp which is obtained 
                         %from independently sample (AA,tAA) pairs
                         
%%                         
for i=10:10
    %-------------------------------------------------------------------------
    %Sample calories, protein, carb and fats in a specific diet
    %-------------------------------------------------------------------------
    Diet=DietList{i};
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
        prob_samp_diet.c=rand(n_foods,1)-0.5;
        [~,res1]=mosekopt('minimize echo(0)',prob_samp_diet);
        [~,res2]=mosekopt('maximize echo(0)',prob_samp_diet);
        x0=[x0 (res1.sol.bas.xx+res2.sol.bas.xx)/2];       
    end
    x0=mean(x0')';
    clear res1 res2
    list_carbfat_samp{i}=ACHR_Sampler(prob_samp_diet,K_CarbFat_all,K_CarbFat_all*x0,nWarmup,nFinal);
end

% Process the sampled values
list_carb_ratio=cell(1,3);
list_lipid_ratio=cell(1,3);
for i=1:3
    list_carb_ratio{i}=zeros(50000,10);
    list_lipid_ratio{i}=zeros(50000,10);
end
for i=1:10
    a=list_carbfat_samp{i};
    x=[a(:,2:3) a(:,1)-a(:,2)-a(:,3)];
    x(x<0)=0;
    x=x./sum(x,2);
    for j=1:3
        b=list_carb_ratio{j};
        b(:,i)=x(:,j);
        list_carb_ratio{j}=b;
    end
    x=a(:,4:6);
    x=x./sum(x,2);
    for j=1:3
        b=list_lipid_ratio{j};
        b(:,i)=x(:,j);
        list_lipid_ratio{j}=b;
    end
end

% Compute and plot average carb/fat composition in 10 pre-defined diets
carb_ratio_diet_average=zeros(10,3);
lipid_ratio_diet_average=zeros(10,3);
for i=1:10
    x=list_carbfat_samp{i};
    carb_ratio_diet_average(i,:) = nanmean([x(:,2:3) x(:,1)-x(:,2)-x(:,3)]./x(:,1));%./sum(mean(x(:,1:3)));
    lipid_ratio_diet_average(i,:) = nanmean(x(:,4:6)./sum(x(:,4:6),2));%./sum(mean(x(:,4:6)));
end

%Plot distribution of carb/fat composition in diets
CarbNames={'Dietary fiber','Sugar','Other'};
FatNames={'Saturated','Polyunsaturated','Monounsaturated'};
figure;
for i=1:3
    x=list_carb_ratio{i};
    subplot(2,3,i);
    [~,idx]=sort(median(x));
    violinplot(x(:,idx),DietNames(idx),'ShowData',false);
    xtickangle(45);
    title(CarbNames{i});
    xlim([0 11]);
    box on;
    
    x=list_lipid_ratio{i};
    subplot(2,3,i+3);
    [~,idx]=sort(median(x));
    violinplot(x(:,idx),DietNames(idx),'ShowData',false);
    xtickangle(45);
    title(FatNames{i});
    xlim([0 11]);
    box on;
end

% Compare ANOVA F-statistic between amino acids and carbs/fats
F_ANOVA_AA_Diet=zeros(18,1);
for i=1:18
    aa_ratio_samp=zeros(nFinal,10);
    for j=1:10
        a=list_aa_samp{j};
        aa_ratio_samp(:,j)=a(:,i)./sum(a,2);
    end
    [~,tbl]=anova1(aa_ratio_samp,[],'off');
    F_ANOVA_AA_Diet(i)=tbl{2,5};   
end

F_ANOVA_Carb_Diet=zeros(3,1);
F_ANOVA_Fat_Diet=zeros(3,1);
for i=1:3
    [~,tbl]=anova1(list_carb_ratio{i},[],'off');
    F_ANOVA_Carb_Diet(i)=tbl{2,5};
    [~,tbl]=anova1(list_lipid_ratio{i},[],'off');
    F_ANOVA_Fat_Diet(i)=tbl{2,5};
end