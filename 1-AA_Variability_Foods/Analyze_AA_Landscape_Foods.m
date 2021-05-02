%% Load food nutrient composition data from the USDA SR database
x=readtable('FoodMatrix.csv','PreserveVariableNames',true);
FoodNames=table2array(x(:,1)); %Name of foods
FoodTypes=table2array(x(:,2)); %Type of foods
FoodMatrix=table2array(x(:,3:end)); %Nutrient values of foods
Nutrients=x.Properties.VariableNames(3:end); %Name of nutrients
load MinServing.txt % Serving size in grams for each food
clear x

%% Create variables related to indice of important nutrients (i.e. nutrients
% considered in USDA nutritional goals and amino acids
load Bounds.txt % Bounds defined by USDA nutritional goals
BdPos=find(max(Bounds')>0); %Indice of nutrients considered in USDA nutritional goals
ProteinPos=140;
CarbohydratePos=34;
SugarPos=82;
FatPos=43;
SatFatPos=157;
AAPos=[5 7 19 22 42 49 50 52 55 81 110 112 113 137 142 168 172 174]; %Indice of amino acids in FoodMatrix

%% Extract foods without missing values in the important nutrients defined in the previous section
CompletePos=find(min(FoodMatrix(:,[BdPos(:);ProteinPos;CarbohydratePos;FatPos;SugarPos;SatFatPos;AAPos']),[],2)~=-1);

%% Characterize variability in AA composition of foods
AANames={'Serine','Tyrosine','Glycine','Phenylalanine','Proline','Valine',...
    'Lysine','Aspartate+Asparagine','Leucine','Isoleucine','Tryptophan','Arginine',...
    'Glutamate+Glutamine','Methionine','Histidine','Threonine','Alanine','Cystine'};

% Select foods with protein (protein weight/dry weight > 0.1)
gProt_gDW=FoodMatrix(CompletePos,ProteinPos)./(100-FoodMatrix(CompletePos,165)); %Fraction of protein in dry weight of food
Foods_wProt=find(gProt_gDW>0.1); %Foods with at least 10%w/dw protein were kept for the following analysis
gAA_gProt=FoodMatrix(CompletePos,AAPos)./sum(FoodMatrix(CompletePos,AAPos),2);
[~,idx]=sort(median(gAA_gProt(Foods_wProt,:)));
figure;violinplot(gAA_gProt(Foods_wProt,idx),AANames(idx),'ShowData',false);
xtickangle(45);
ylabel('Amino acid [g/gProtein]');
title('Distributions of amino acid composition in foods');

%% Remove categories with too few protein-rich foods
CATs=FoodTypes(CompletePos(Foods_wProt));
CATs(strcmp(CATs,'Meals, Entrees, and Side Dishes')...
    | strcmp(CATs,'Soups, Sauces, and Gravies')...
    | strcmp(CATs,'Restaurant Foods'))={'Others'}; 
CATs(strcmp(CATs,'Nut and Seed Products'))={'Nuts'};
CATs(strcmp(CATs,'Breakfast Cereals')...
    | strcmp(CATs,'Cereal Grains and Pasta'))={'Cereals'};
CATs(strcmp(CATs,'Sausages and Luncheon Meats'))={'Others'};
CATs(strcmp(CATs,'Baked Goods'))={'Others'};
unqCATs=unique(CATs);
for i=1:length(unqCATs)
    n=sum(strcmp(unqCATs{i},CATs));
    if n<20
        CATs(strcmp(unqCATs{i},CATs)==true)={'Others'};
    end
end
map0 = brewermap(12,'Set3');
map1 = brewermap(3,'Accent');
map = [map0;map1];
clear map0 map1

%% tSNE analysis of AA profiles in foods
x=tsne(gAA_gProt(Foods_wProt,:));
figure;gscatter(x(:,1),x(:,2),CATs,map,[],10);xlabel('tSNE 1');ylabel('tSNE 2');
title('t-SNE analysis of amino acid profiles in foods');
clear x

%% PCA of amino acid profiles in foods
[~,score,~,~,explained,~]=pca(gAA_gProt(Foods_wProt,:));
figure;gscatter(score(:,1),score(:,2),CATs,map,[],10);
title('PCA of amino acid profiles in foods');

%% Coefficient of variation
CV_AA=std(gAA_gProt(Foods_wProt,:))./mean(gAA_gProt(Foods_wProt,:));

%% Average amino acid content in each category
CATs_keep=setdiff(unique(CATs),'Others');
MeanAA_CAT=zeros(length(CATs_keep),18);
CVAA_CAT=zeros(length(CATs_keep),18);
%CATs=FoodTypes(CompletePos(Foods_wProt));
CAT_Count=zeros(length(CATs_keep),1);
for i=1:length(CATs_keep)
    idx=find(strcmp(CATs_keep{i},CATs)==true);
    CAT_Count(i)=length(idx);
    if length(idx)>1
        MeanAA_CAT(i,:)=mean(gAA_gProt(Foods_wProt(idx),:));
        CVAA_CAT(i,:)=std(gAA_gProt(Foods_wProt(idx),:))./MeanAA_CAT(i,:);
    else
        MeanAA_CAT(i,:)=gAA_gProt(Foods_wProt(idx),:);
        CVAA_CAT(i,:)=0;
    end
end
figure;plotFractions(MeanAA_CAT,CATs_keep,AANames);
title('Amino acid composition of foods');
figure;heatmap_cluster(CVAA_CAT,CATs_keep,AANames,[0 0.8]);
title('CV of amino acid composition in food categories');

% Plot heatmap to compare average AA abundance in food groups
data = (MeanAA_CAT-min(MeanAA_CAT))./(max(MeanAA_CAT)-min(MeanAA_CAT));
figure;
cmap_now = brewermap(100,'RdYlGn');
colormap(cmap_now(75:-1:26,:));
heatmap_cluster(data',AANames,CATs_keep,[0 1]);
clear cmap_now data
colorbar('Ticks',[0 1],'TickLabels',{'0','Row max'});
title('Amino acid abundances in food groups');

%% Compare variability in AA composition with variability of fats and carbs

% Extract carbohydrate/fat content information in foods
CarbNames={'Carbohydrate, by difference,g','Fiber, total dietary,g',...
    'Sugars, total,g'};
[~,pos]=ismember(CarbNames,Nutrients);
CarbNames={'Dietary fiber','Sugar','Other'};
x=FoodMatrix(:,pos);x(x==-1)=NaN;
CarbRatio=[x(:,2:3) x(:,1)-x(:,2)-x(:,3)];
CarbRatio(CarbRatio<0)=0;
CarbRatio=CarbRatio./sum(CarbRatio,2);

LipidNames={'Total lipid (fat),g','Fatty acids, total saturated,g',...
'Fatty acids, total polyunsaturated,g','Fatty acids, total monounsaturated,g',...
'18:1 undifferentiated,g','18:2 undifferentiated,g',...
'18:3 undifferentiated,g','16:0,g','18:0,g'};
[~,pos]=ismember(LipidNames,Nutrients);
LipidNames_1={'Saturated','Polyunsaturated','Monounsaturated'};
LipidNames_2={'18:1','18:2','18:3','16:0','18:0','Other'};
LipidNames_Combine = [LipidNames_1 LipidNames_2];
x=FoodMatrix(:,pos);x(x==-1)=NaN;
LipidRatio_1=x(:,2:4);LipidRatio_1=LipidRatio_1./x(:,1);
LipidRatio_2=[x(:,5:9) x(:,1)-sum(x(:,5:9),2)];
LipidRatio_2(LipidRatio_2<0)=0;LipidRatio_2=LipidRatio_2./x(:,1);
LipidRatio_Combine = [LipidRatio_1 LipidRatio_2];

% One-way ANOVA comparing AA/carb/fat composition across foods
FoodTypes_Reduced=FoodTypes(CompletePos);
idx_keep = find(ismember(FoodTypes_Reduced,CATs_keep));

F_ANOVA_AA=zeros(18,1);
F_ANOVA_Carb=zeros(3,1);
F_ANOVA_Fat=zeros(9,1);
for i=1:18
    [~,tbl]=anova1(gAA_gProt(idx_keep,i),FoodTypes_Reduced(idx_keep),'off');
    F_ANOVA_AA(i)=tbl{2,5};
end
for i=1:3
    [~,tbl]=anova1(CarbRatio(CompletePos(idx_keep),i),FoodTypes_Reduced(idx_keep),'off');
    F_ANOVA_Carb(i)=tbl{2,5};
end
for i=1:9
    [~,tbl]=anova1(LipidRatio_Combine(CompletePos(idx_keep),i),FoodTypes_Reduced(idx_keep),'off');
    F_ANOVA_Fat(i)=tbl{2,5};
end

% Show violin plots comparing 18:0, polyunsaturated, dietary fiber, lysine,
% methionine, histidine across food groups
data_violin = [LipidRatio_Combine(CompletePos(idx_keep),[2 8 1]) CarbRatio(CompletePos(idx_keep),1) gAA_gProt(idx_keep,[14 15 7 5])]; 
title_violin = {'Polyunsaturated fat','Saturated fat 18:0','Saturated fat',...
    'Dietary fiber','Methionine','Histidine','Lysine','Proline'};
ylabel_violin = {'g/g total fat','g/g total fat','g/g total fat','g/g total carbohydrate',...
    'g/g total amino acids','g/g total amino acids','g/g total amino acids','g/g total amino acids'};
data_violin(isnan(data_violin)) = 0;
figure;
for i=1:8
    subplot(2,4,i);
    violinplot(data_violin(:,i),FoodTypes_Reduced(idx_keep),'ViolinColor',[112 173 71]/255);
    title(title_violin(i));
    xtickangle(45);
    ylabel(ylabel_violin(i));
    box on;
end

clear data_violin title_violin ylabel_violin
