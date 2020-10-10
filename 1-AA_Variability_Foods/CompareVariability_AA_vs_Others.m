%Compare variability in AA composition with variability of fats and carbs

%% Extract carbohydrate/fat content information in foods
CarbNames={'"Carbohydrate, by difference",g','"Fiber, total dietary",g',...
    '"Sugars, total",g'};
[~,pos]=ismember(CarbNames,Nutrients);
CarbNames={'Dietary fiber','Sugar','Other'};
x=FoodMatrix(:,pos);x(x==-1)=NaN;
CarbRatio=[x(:,2:3) x(:,1)-x(:,2)-x(:,3)];
CarbRatio(CarbRatio<0)=0;
CarbRatio=CarbRatio./sum(CarbRatio,2);

LipidNames={'"Total lipid (fat)",g','"Fatty acids, total saturated",g',...
'"Fatty acids, total polyunsaturated",g','"Fatty acids, total monounsaturated",g',...
'"18:1 undifferentiated",g','"18:2 undifferentiated",g',...
'"18:3 undifferentiated",g','"16:0",g','"18:0",g'};
[~,pos]=ismember(LipidNames,Nutrients);
LipidNames_1={'Saturated','Polyunsaturated','Monounsaturated'};
LipidNames_2={'18:1','18:2','18:3','16:0','18:0','Other'};
LipidNames_Combine = [LipidNames_1 LipidNames_2];
x=FoodMatrix(:,pos);x(x==-1)=NaN;
LipidRatio_1=x(:,2:4);LipidRatio_1=LipidRatio_1./sum(LipidRatio_1,2);
LipidRatio_2=[x(:,5:9) x(:,1)-sum(x(:,5:9),2)];
LipidRatio_2(LipidRatio_2<0)=0;LipidRatio_2=LipidRatio_2./sum(LipidRatio_2,2);
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

%%
% Show violin plots comparing 18:0, polyunsaturated, dietary fiber, lysine,
% methionine, histidine across food groups
data_violin = [LipidRatio_Combine(CompletePos(idx_keep),[8 2]) CarbRatio(CompletePos(idx_keep),1) gAA_gProt(idx_keep,[15 14 7])]; 
title_violin = {'18:0','Polyunsaturated','Dietary fiber','Lysine','Methionine','Histidine'};
data_violin(isnan(data_violin)) = 0;
figure;
for i=1:6
    subplot(2,3,i);
    violinplot(data_violin(:,i),FoodTypes_Reduced(idx_keep));
    title(title_violin(i));
    xtickangle(45);
    box on;
end

clear data_violin title_violin

%% Part II: variability in human dietary records
[~,pos]=ismember({'Carbohydrate(gm)','Dietary fiber(gm)','Total sugars(gm)'},NutAlias);
x=Mat_Nutrients(:,pos);
CarbRatio_NHANES=[x(:,2:3) x(:,1)-x(:,2)-x(:,3)];
CarbRatio_NHANES(CarbRatio_NHANES<0)=0;
CarbRatio_NHANES=CarbRatio_NHANES./sum(CarbRatio_NHANES,2);

[~,pos]=ismember({'Total saturated fatty acids(gm)','Total polyunsaturated fatty acids(gm)',...
    'Total monounsaturated fatty acids(gm)'},NutAlias);
x=Mat_Nutrients(:,pos);
LipidRatio_NHANES=x./sum(x,2);

%% Part III: difference between processed and unprocessed foods
Raw_Start=FoodNames(contains(FoodNames,' raw'));
Raw_All=union(Raw_Start,Raw_Add);
Processed_All=setdiff(FoodNames,union(Raw_All,Culinary_Ingredients));
[~,pos_raw]=ismember(Raw_All,FoodNames); %Raw foods
[~,pos_processed]=ismember(Processed_All,FoodNames); %Processed foods
[~,pos_ci]=ismember(Culinary_Ingredients,FoodNames); %Culinary Ingredients (oils, spices, herbs etc)
ProcessTag=zeros(8788,1);
ProcessTag(pos_ci)=1;
ProcessTag(pos_processed)=2;

figure;
for i=1:18
    subplot(3,6,i);
    violinplot(gAA_gProt(ProcessTag(CompletePos)~=1,i),...
        ProcessTag(CompletePos(ProcessTag(CompletePos)~=1)),...
        'ShowData',false);
    title(AANames(i));
    xticklabels({'Raw','Processed'});xtickangle(45);box on;
    pAA_RawProcessed(i)=ranksum(gAA_gProt(ProcessTag(CompletePos)==0,i),...
        gAA_gProt(ProcessTag(CompletePos)==2,i));
end
figure;
for i=1:3
    subplot(2,3,i);
    violinplot(CarbRatio(ProcessTag~=1,i),...
        ProcessTag(ProcessTag~=1),'ShowData',false);
    title(CarbNames(i));
    xticklabels({'Raw','Processed'});xtickangle(45);box on;
    pCarb_RawProcessed(i)=ranksum(CarbRatio(ProcessTag==0,i),...
        CarbRatio(ProcessTag==2,i));
end
for i=1:3
    subplot(2,3,i+3);
    violinplot(LipidRatio_1(ProcessTag~=1,i),...
        ProcessTag(ProcessTag~=1),'ShowData',false);
    title(LipidNames_1(i));
    xticklabels({'Raw','Processed'});xtickangle(45);box on;
    pFat_RawProcessed(i)=ranksum(LipidRatio_1(ProcessTag==0,i),...
        LipidRatio_1(ProcessTag==2,i));
end

%% Part IV: Compare raw and processed foods for each category
[tbl,chi2,p,labels]=crosstab(ProcessTag(CompletePos),FoodTypes(CompletePos));
FoodCat_RawProcessed=labels(tbl(1,:)>1 & tbl(3,:)>1,2); %Food categories that include both raw and processed foods
pAA_RawProcessed_FoodCat=zeros(length(FoodCat_RawProcessed),18);
logfcAA_RawProcessed_FoodCat=zeros(length(FoodCat_RawProcessed),18);

pCarb_RawProcessed_FoodCat=zeros(length(FoodCat_RawProcessed),3);
logfcCarb_RawProcessed_FoodCat=zeros(length(FoodCat_RawProcessed),3);

pFat_RawProcessed_FoodCat=zeros(length(FoodCat_RawProcessed),3);
logfcFat_RawProcessed_FoodCat=zeros(length(FoodCat_RawProcessed),3);
for i=1:length(FoodCat_RawProcessed)
    pos_raw=find(strcmp(FoodTypes(CompletePos),FoodCat_RawProcessed{i}) & ProcessTag(CompletePos)==0);
    pos_processed=find(strcmp(FoodTypes(CompletePos),FoodCat_RawProcessed{i}) & ProcessTag(CompletePos)==2);
    logfcAA_RawProcessed_FoodCat(i,:)=log(median(gAA_gProt(pos_processed,:))./median(gAA_gProt(pos_raw,:)));
    logfcCarb_RawProcessed_FoodCat(i,:)=log(nanmedian(CarbRatio(CompletePos(pos_processed),:))./nanmedian(CarbRatio(CompletePos(pos_raw),:)));
    logfcFat_RawProcessed_FoodCat(i,:)=log(nanmedian(LipidRatio_1(CompletePos(pos_processed),:))./nanmedian(LipidRatio_1(CompletePos(pos_raw),:)));
    for j=1:18
        pAA_RawProcessed_FoodCat(i,j)=ranksum(gAA_gProt(pos_processed,j),gAA_gProt(pos_raw,j));
        if j<=3
            pCarb_RawProcessed_FoodCat(i,j)=ranksum(CarbRatio(CompletePos(pos_processed),j),...
                CarbRatio(CompletePos(pos_raw),j));
            pFat_RawProcessed_FoodCat(i,j)=ranksum(LipidRatio_1(CompletePos(pos_processed),j),...
                LipidRatio_1(CompletePos(pos_raw),j));
        end
    end
end
logfcAA_RawProcessed_FoodCat(pAA_RawProcessed_FoodCat>0.05)=0;
logfcCarb_RawProcessed_FoodCat(pCarb_RawProcessed_FoodCat>0.05)=0;
logfcCarb_RawProcessed_FoodCat(isnan(logfcCarb_RawProcessed_FoodCat) | isinf(logfcCarb_RawProcessed_FoodCat))=0;
logfcFat_RawProcessed_FoodCat(pFat_RawProcessed_FoodCat>0.05)=0;
figure;
%subplot(1,3,1);
heatmap_cluster(logfcAA_RawProcessed_FoodCat,FoodCat_RawProcessed,AANames,[-0.4 0.4]);
%{
subplot(1,3,2);
heatmap_cluster(logfcCarb_RawProcessed_FoodCat(:,1:2),FoodCat_RawProcessed,CarbNames(1:2),[-1 1]);
subplot(1,3,3);
heatmap_cluster(logfcFat_RawProcessed_FoodCat,FoodCat_RawProcessed,LipidNames_1,[-1 1]);
%}

%% Part V: Variability between pre-defined diets
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