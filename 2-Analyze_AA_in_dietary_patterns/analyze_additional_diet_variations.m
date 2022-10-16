%-------------------------------------------------------------------------
% Define some additional variations of diets
%-------------------------------------------------------------------------
addpath(genpath('../functions'));
load ../data/AA_variables.mat;
define_diets;

%-------------------------------------------------------------------------
%Additional diet 1: Atkins diet at the "ongoing weight loss" phase
%-------------------------------------------------------------------------
%allows less than 50g carbohydrates per day
%-------------------------------------------------------------------------
ConstraintNames={'Carbohydrate'};
ConstraintLBs=0;
ConstraintUBs=50;
ConstraintMatrix=0.01*FoodMatrix(:,34)';
%-------------------------------------------------------------------------
%Construct the data structure
%-------------------------------------------------------------------------
Atkins_V1.ConstraintNames=ConstraintNames;
Atkins_V1.ConstraintLBs=ConstraintLBs;
Atkins_V1.ConstraintUBs=ConstraintUBs;
Atkins_V1.ConstraintMatrix=ConstraintMatrix;
%-------------------------------------------------------------------------
%Finished definition of Atkins diet, variation 1
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
%Additional diet 2: The Eco-Atkins diet which is a combination of Atkins
%and vegan (plant-based) diet
%-------------------------------------------------------------------------
%allows less than 20g carbohydrates per day
%allows no consumption of animal products
%-------------------------------------------------------------------------
ConstraintNames={'Carbohydrate', 'Meat, dairy and eggs'};
ConstraintLBs=[0;0];
ConstraintUBs=[20;0];
ConstraintMatrix(1,:)=0.01*FoodMatrix(:,34)';
ConstraintMatrix(2,:)=double(ismember(FoodTypeIDs',[1 2 3 4 8 9 10 13 15 17 18 19 20 21 22 23]));
animal_fats = readtable('input_data/list_animal_dairy_fat.csv'); % Read the list of animal and dairy products in the category "fat"
ConstraintMatrix(2,:) = ConstraintMatrix(2,:) + double(ismember(FoodNames,animal_fats.Name))';
%-------------------------------------------------------------------------
%Construct the data structure
%-------------------------------------------------------------------------
Atkins_V2.ConstraintNames=ConstraintNames;
Atkins_V2.ConstraintLBs=ConstraintLBs;
Atkins_V2.ConstraintUBs=ConstraintUBs;
Atkins_V2.ConstraintMatrix=ConstraintMatrix;
%-------------------------------------------------------------------------
%Finished definition of Atkins diet, variation 2
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
%Additional diet 3: The Atkins diet with limited consumption of red meat
%-------------------------------------------------------------------------
%allows less than 20g carbohydrates per day
%allows no consumption of red meat
%-------------------------------------------------------------------------
ConstraintNames={'Carbohydrate', 'Red meat'};
ConstraintLBs = [0;0];
ConstraintUBs = [20;0];
ConstraintMatrix(1,:) = 0.01*FoodMatrix(:,34)';
ConstraintMatrix(2,:) = double(ismember(FoodTypeIDs',[4 13 17 20]));
%-------------------------------------------------------------------------
%Construct the data structure
%-------------------------------------------------------------------------
Atkins_V3.ConstraintNames=ConstraintNames;
Atkins_V3.ConstraintLBs=ConstraintLBs;
Atkins_V3.ConstraintUBs=ConstraintUBs;
Atkins_V3.ConstraintMatrix=ConstraintMatrix;
%-------------------------------------------------------------------------
%Finished definition of Atkins diet, variation 3
%-------------------------------------------------------------------------

%-----------------------------------------------------------------------
%Additional Diet 4: The traditional Mediterranean diet
%Recommended composition of a Mediterranean diet is taken from the
%Dietary guidelines for adults in Greece. 
%Source: https://www.mednet.gr/archives/1999-5/pdf/516.pdf
%-----------------------------------------------------------------------
ConstraintNames={'Whole-grain cereals','Fish','Fruits','Poultry',...
    'Sweets','Vegetables','Potatoes','Red meat',...
    'Olives/nuts/seeds','Olive oil','Eggs','Dairy','Others'};
% Read list of refined cereals
refined_cereals = readtable('input_data/list_refined_cereals.csv');


%Lower and upper bounds of servings of foods in each category
ConstraintLBs=[8;5/7;3;4/7;0;6;3/7;0;3/7;1;3/7;2;0];
ConstraintUBs=[8;6/7;3;4/7;3/7;6;3/7;4/30;4/7;Inf;3/7;2;0];
%Coefficient matrix for the constraints
ConstraintMatrix=zeros(13,8788);
%Cereals
ConstraintMatrix(1,:)=ReciprocalServing'.*(FoodTypeIDs'==7);
ConstraintMatrix(1,ismember(FoodNames,refined_cereals.Name)) = 0; %Remove refined cereals
%Fish
ConstraintMatrix(2,:)=ReciprocalServing'.*(FoodTypeIDs'==21);
%Fruits
ConstraintMatrix(3,:)=ReciprocalServing'.*(FoodTypeIDs'==12 & contains(FoodNames','juice')==0);
%Poultry
ConstraintMatrix(4,:)=ReciprocalServing'.*(FoodTypeIDs'==18);
%Sweets
ConstraintMatrix(5,:)=ReciprocalServing'.*(FoodTypeIDs'==25);
%Vegetables or legumes
ConstraintMatrix(6,:)=ReciprocalServing'.*(FoodTypeIDs'==26 | FoodTypeIDs' == 14);
%Potatoes
ConstraintMatrix(7,:)=ReciprocalServing'.*(contains(FoodNames','Potatoes'));
%Red meat
ConstraintMatrix(8,:)=ReciprocalServing'.*(FoodTypeIDs'==4 | FoodTypeIDs'==13 | FoodTypeIDs'==17);
%Olives/nuts/seeds
ConstraintMatrix(9,:)=ReciprocalServing'.*(FoodTypeIDs'==16 | contains(FoodNames','Olives'));
%Olive oil
ConstraintMatrix(10,:)=ReciprocalServing'.*(contains(FoodNames','Oil, olive'));
%Eggs
ConstraintMatrix(11,:)=ReciprocalServing'.*(FoodTypeIDs'==9);
%Dairy
ConstraintMatrix(12,:)=ReciprocalServing'.*(FoodTypeIDs'==8);
%Others
ConstraintMatrix(13,:)=(max(ConstraintMatrix(1:12,:))==0);
%-------------------------------------------------------------------------
%Construct the data structure
%-------------------------------------------------------------------------
Mediterranean_V1.ConstraintNames=ConstraintNames;
Mediterranean_V1.ConstraintLBs=ConstraintLBs;
Mediterranean_V1.ConstraintUBs=ConstraintUBs;
Mediterranean_V1.ConstraintMatrix=ConstraintMatrix;
%-------------------------------------------------------------------------
%Finished definition of Mediterranean diet
%-------------------------------------------------------------------------


%-----------------------------------------------------------------------
%Additional Diet 5: Diet with highest Mediterranean diet score
%The individual should have score 1 on all nine components of the
%Mediterranean diet score (MDS) defined in the paper "Adherence to a
%Mediterranean diet and survival in a greek population"
%-----------------------------------------------------------------------
ConstraintNames={'Vegetables', 'Legumes', 'Fruits and Nuts',...
    'Cereals', 'Fish', 'Meat', 'Dairy', 'Alcohol',...
    'Monounsaturated-to-saturated Ratio'};

%Lower and upper bounds of servings of foods in each category
ConstraintLBs=[499.6;6.7;356.3;139.7;18.8;0;0;5;0];
ConstraintUBs=[Inf;Inf;Inf;Inf;Inf;89.8;191.1;25;Inf];
%Coefficient matrix for the constraints
ConstraintMatrix=zeros(9,8788);
%Vegetables
ConstraintMatrix(1,:)=double(FoodTypeIDs'==26);
%Legumes
ConstraintMatrix(2,:)=double(FoodTypeIDs' == 14);
%Fruits and Nuts
ConstraintMatrix(3,:)=double((FoodTypeIDs'==12 & contains(FoodNames','juice')==0)...
    | (FoodTypeIDs'==16));
%Cereals
ConstraintMatrix(4,:)=double(FoodTypeIDs'==7);
%Fish
ConstraintMatrix(5,:)=double(FoodTypeIDs'==21);
%Meat
ConstraintMatrix(6,:)=double(ismember(FoodTypeIDs',[4 13 17 18 20]));
%Dairy
ConstraintMatrix(7,:)=double(FoodTypeIDs'==8);
%Alcohol
ConstraintMatrix(8,:) = 0.01*FoodMatrix(:,59)';
%Monounsaturated-to-saturated fat ratio
ConstraintMatrix(9,:) = FoodMatrix(:,14)' - 1.7*FoodMatrix(:,157)';

%-------------------------------------------------------------------------
%Construct the data structure
%-------------------------------------------------------------------------
Mediterranean_V2.ConstraintNames=ConstraintNames;
Mediterranean_V2.ConstraintLBs=ConstraintLBs;
Mediterranean_V2.ConstraintUBs=ConstraintUBs;
Mediterranean_V2.ConstraintMatrix=ConstraintMatrix;
%-------------------------------------------------------------------------
%Finished definition of Mediterranean diet
%-------------------------------------------------------------------------

% Make the list of the diet variations
Diet_Variations = {Atkins, Atkins_V1, Atkins_V2, Atkins_V3,...
    Mediterranean, Mediterranean_V1, Mediterranean_V2};

Diet_Variation_Names = {'Atkins,original', 'Atkins,ongoing weight loss',...
    'Eco-Atkins', 'Atkins, restricted red meat',...
    'Mediterranean,pyramid-based', 'Mediterranean, traditional',...
    'Mediterranean,MDS-based'};


%% Calculate ranges of AA daily intake in different diets with calories
% level [1800,2200]

VarDiet_AAmin=zeros(7,18);
VarDiet_AAmax=zeros(7,18);
label_baby_food = double(contains(FoodNames(CompletePos),'Babyfood'))';
for i=1:7
    Diet=Diet_Variations{i};
    Diet_Variation_Names{i}
    prob_aa_diet.blx=zeros(2335,1);
    prob_aa_diet.bux=1000*ones(2335,1); 
    %No more than 1kg for one food
    prob_aa_diet.blc=[Diet.ConstraintLBs;0;-Inf]; 
    prob_aa_diet.buc=[Diet.ConstraintUBs;3000;0];  
    %No more than 3kg total food
    prob_aa_diet.a=[Diet.ConstraintMatrix(:,CompletePos);ones(1,length(CompletePos));label_baby_food];
    %No baby food
    
    % 1800-2200 kcal energy per day
    prob_aa_diet.blc=[prob_aa_diet.blc;1800];
    prob_aa_diet.buc=[prob_aa_diet.buc;2200];
    prob_aa_diet.a=[prob_aa_diet.a;FoodMatrix(CompletePos,143)'/100];
    
    for j=1:18
        c=FoodMatrix(CompletePos,AAPos(j))/100;
        prob_aa_diet.c=c;
        [~,res_min] = mosekopt('minimize echo(0)',prob_aa_diet);
        VarDiet_AAmin(i,j)=res_min.sol.bas.pobjval;
        [~,res_max] = mosekopt('maximize echo(0)',prob_aa_diet);
        VarDiet_AAmax(i,j)=res_max.sol.bas.pobjval;
    end
end

VarDiet_AAmin_Norm=VarDiet_AAmin./max(VarDiet_AAmin); %Normalize to column maximum
VarDiet_AAmax_Norm=VarDiet_AAmax./max(VarDiet_AAmax);
VarDiet_AAVar=log10(VarDiet_AAmax./VarDiet_AAmin); %Variability of daily AA intake
VarDiet_AAVar(isinf(VarDiet_AAVar))=0;
VarDiet_Var_Pos=find(max(VarDiet_AAVar')>0);
VarDiet_AAVar=VarDiet_AAVar(max(VarDiet_AAVar')>0,:);

%Plot ranges of absolute AA intake in human dietary patterns
cmap_now = brewermap(100,'RdYlGn');
figure;
colormap(cmap_now(75:-1:26,:));
subplot(1,2,1);
heatmap_cluster(VarDiet_AAmin_Norm',AANames,Diet_Variation_Names,[0 1]);
title('Minimal daily intake');
subplot(1,2,2);
heatmap_cluster(VarDiet_AAmax_Norm',AANames,Diet_Variation_Names,[0 1]);
title('Maximal daily intake');
colorbar('Ticks',[0 1],'TickLabels',{'0','Row max'});

figure;
colormap(cmap_now(75:-1:26,:));
heatmap_cluster(VarDiet_AAVar',AANames,Diet_Variation_Names(VarDiet_Var_Pos),[0 3]);
title('Variability');
colorbar('Ticks',[0 3]);

clear cmap_now

%% Sample random diets under the additional variations of Med and Atkins diet

% Define variables
%Coefficients for sampling AA compositions
K_AA_all=FoodMatrix(CompletePos,AAPos)'/100;
%K_totalAA=sum(K_AA_all);
nWarmup=1000;
nFinal=50000;

% Sample AA composition of diets
list_aa_samp_var_diets=cell(1,7); %Used to be list_aa_ratio_samp which is obtained 
                         %from independently sample (AA,tAA) pairs
for i=1:7
    Diet=Diet_Variations{i};
    Diet_Variation_Names(i)
    prob_samp_diet.blx=zeros(2335,1);
    prob_samp_diet.bux=1000*ones(2335,1); %No more than 1kg for one food
    prob_samp_diet.blc=[Diet.ConstraintLBs;0;-Inf];
    prob_samp_diet.buc=[Diet.ConstraintUBs;3000;0];  %No more than 3kg total food
    prob_samp_diet.a=[Diet.ConstraintMatrix(:,CompletePos);ones(1,length(CompletePos));label_baby_food];
    %No baby food
    prob_samp_diet.blc=[prob_samp_diet.blc;1800];
    prob_samp_diet.buc=[prob_samp_diet.buc;2200];
    prob_samp_diet.a=[prob_samp_diet.a;FoodMatrix(CompletePos,143)'/100];
    
    x0=[];
    for j=1:100
        prob_samp_diet.c=rand(2335,1)-0.5;
        [~,res1]=mosekopt('minimize echo(0)',prob_samp_diet);
        [~,res2]=mosekopt('maximize echo(0)',prob_samp_diet);
        x0=[x0 (res1.sol.bas.xx+res2.sol.bas.xx)/2];       
    end
    x0=mean(x0')';
    clear res1 res2
    list_aa_samp_var_diets{i}=ACHR_Sampler(prob_samp_diet,K_AA_all,K_AA_all*x0,nWarmup,nFinal);
end

%% Process the sampled values and show violin plots

% Index for Atkins and Mediterranean diets
Atkins_idx = [1 2 4];
Mediterranean_idx = [5 6 7];
F_ANOVA_Diet_Vars = zeros(18,2);

fh1 = figure;
fh2 = figure;
for i=1:18
    aa_ratio_samp=zeros(nFinal,7);
    for j=1:7
        a=list_aa_samp_var_diets{j};
        aa_ratio_samp(:,j)=a(:,i)./sum(a,2);
    end
    %Compute ANOVA F-statistic
    [~,tbl] = anova1(aa_ratio_samp(:,Atkins_idx),1:3,'off');
    F_ANOVA_Diet_Vars(i,1) = tbl{2,5};
    [~,tbl] = anova1(aa_ratio_samp(:,Mediterranean_idx),1:3,'off');
    F_ANOVA_Diet_Vars(i,2) = tbl{2,5};
    
    figure(fh1);
    subplot(3,6,i);
    violinplot(aa_ratio_samp(:,Atkins_idx),...
        Diet_Variation_Names(Atkins_idx),'ShowData',false);
    xtickangle(45);
    if i < 13
        xticklabels([]);
    end
    title(AANames{i});
    xlim([0 4]);
    box on;
    
    figure(fh2);
    subplot(3,6,i);
    violinplot(aa_ratio_samp(:,Mediterranean_idx),...
        Diet_Variation_Names(Mediterranean_idx),'ShowData',false);
    xtickangle(45);
    if i < 13
        xticklabels([]);
    end
    title(AANames{i});
    xlim([0 4]);
    box on;
end
mean_aa_ratio_var_diets=zeros(7,18);
for i=1:7
    a=list_aa_samp_var_diets{i};
    mean_aa_ratio_var_diets(i,:)=mean(a./sum(a,2));
end

% Plot average fractions of AAs in dietary patterns
data = (mean_aa_ratio_var_diets-min(mean_aa_ratio_var_diets))./(max(mean_aa_ratio_var_diets)-min(mean_aa_ratio_var_diets));
figure;

cmap_now = brewermap(100,'RdYlGn');
heatmap_cluster(data',AANames,Diet_Variation_Names,[0 1],cmap_now(75:-1:26,:));
colorbar('Ticks',[0 1],'TickLabels',{'Row min','Row max'});
title('Relative abundance of amino acids in diets');
clear data cmap_now

% PCA of amino acid compositions of diets
diet_tag=reshape(repmat(Diet_Variation_Names,5000,1),35000,1);
aa_ratio_downsample_var_diets=zeros(35000,18);
for i=1:7
    rp=randperm(35000);
    a=list_aa_samp_var_diets{i};
    aa_ratio_downsample_var_diets((i-1)*5000+1:i*5000,:)=a(rp(1:5000),:)./sum(a(rp(1:5000),:),2);    
end

[~,score,~,~,explained,~]=pca(aa_ratio_downsample_var_diets);
figure;
color=brewermap(11,'Set3');
for i=1:7
    data=score((i-1)*5000+1:i*5000,:);
    scatter3(data(:,1),data(:,2),data(:,3),5,color(i,:),'filled');
    hold(gca,'on');
end
box on;
grid on;
xlabel(sprintf('PC1 (%.1f%%)',explained(1)));
ylabel(sprintf('PC2 (%.1f%%)',explained(2)));
zlabel(sprintf('PC3 (%.1f%%)',explained(3)));
legend(Diet_Variation_Names);
title('PCA of relative amino acid compositions of diets');

%% Compare additional variations of diets with the original ones in terms of mean AAs
mean_aa_ratio_comb = [mean_aa_ratio;mean_aa_ratio_var_diets([2 4 6 7],:)];
diet_names_comb = [DietNames Diet_Variation_Names([2 4 6 7])];
data = (mean_aa_ratio_comb-min(mean_aa_ratio_comb))./(max(mean_aa_ratio_comb)-min(mean_aa_ratio_comb));
figure;

cmap_now = brewermap(100,'RdYlGn');
heatmap_cluster(data',AANames,diet_names_comb,[0 1],cmap_now(75:-1:26,:));
colorbar('Ticks',[0 1],'TickLabels',{'Row min','Row max'});
title('Relative abundance of amino acids in diets');
clear data cmap_now