%% Load food nutrient composition data from the USDA SR database
x=readtable('FoodMatrix.csv');
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

%% Generate heatmap as output
data = (MeanAA_CAT-min(MeanAA_CAT))./(max(MeanAA_CAT)-min(MeanAA_CAT));
figure;
heatmap_cluster(data',AANames,CATs_keep,[0 1]);
cmap_now = brewermap(100,'PiYG');
colormap(cmap_now(75:-1:26,:));
clear cmap_now data
colorbar('Ticks',[0 1],'TickLabels',{'0','Row max'});
title('Amino acid abundances in food groups');
