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
x=FoodMatrix(:,pos);x(x==-1)=NaN;
LipidRatio_1=x(:,2:4);LipidRatio_1=LipidRatio_1./sum(LipidRatio_1,2);
LipidRatio_2=[x(:,5:9) x(:,1)-sum(x(:,5:9),2)];
LipidRatio_2(LipidRatio_2<0)=0;LipidRatio_2=LipidRatio_2./sum(LipidRatio_2,2);



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