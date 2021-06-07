%--------------------------------------------------------------------------
% Analyze the association between obesity and dietary amino acid
% composition in detail
%--------------------------------------------------------------------------

% Add necessary functions
addpath(genpath('../functions'));

%% Load existing variables from workspace saved in ../data
load ../data/AA_variables.mat;

%% Plot heatmap for association between dietary amino acid composition and obesity
NegAANames={'Phenylalanine','Aspartate+Asparagine','Tryptophan','Valine',...
    'Glutamate+Glutamine'};
PosAANames={'Glycine','Alanine','Methionine','Lysine','Histidine'};
UshapeAANames={'Threonine','Proline','Arginine','Leucine','Cystine'};
[~,NegAAPos_NHANES]=ismember(NegAANames,AANames_NHANES);
[~,PosAAPos_NHANES]=ismember(PosAANames,AANames_NHANES);
[~,UAAPos_NHANES]=ismember(UshapeAANames,AANames_NHANES);
nq=5;
ratios=zeros(18,nq);
label_obesity=DiseaseScores(:,3);
goodpos=find(~isnan(label_obesity)); %Indice for individuals with BMI record
label_obesity=label_obesity(goodpos);
p_Chi2_AA_Obesity=zeros(1,18);
for i=1:18
    x=Var_Diet_Adjusted(goodpos,39+i);
    q4=quantile(x,0:1/nq:1);
    label_quantile=zeros(length(goodpos),1);
    for k=1:nq
        qpos=find(x>=q4(k) & x<q4(k+1));
        if k==nq
            qpos=find(x>=q4(k));
        end
        label_quantile(qpos)=k;
        ratios(i,k)=sum(DiseaseScores(goodpos(qpos),3)==1)/length(qpos);
    end
    [~,~,p_Chi2_AA_Obesity(i)]=crosstab(label_obesity,label_quantile);
end

figure;
subplot(3,1,1);
heatmap(ratios(NegAAPos_NHANES,:),1:5,NegAANames,[],'MinColorValue',0.34,...
    'MaxColorValue',0.44,'ShowAllTicks',true,'GridLine',':');
title('Negative association');
subplot(3,1,2);
heatmap(ratios(PosAAPos_NHANES,:),1:5,PosAANames,[],'MinColorValue',0.34,...
    'MaxColorValue',0.44,'ShowAllTicks',true,'GridLine',':');
title('Positive association');
subplot(3,1,3);
heatmap(ratios(UAAPos_NHANES,:),1:5,UshapeAANames,[],'MinColorValue',0.34,...
    'MaxColorValue',0.44,'ShowAllTicks',true,'GridLine',':');
title('U-shaped relationship');
xlabel('Dietary intake quantile');
map=brewermap(64,'RdBu');map=map(64:-1:1,:);colormap(map);colorbar;


%% Plot combinational features consisting of AAs with positive or negative
%association with obesity
AACodes_NHANES={'W','T','I','L','K','M','C','F','Y','V','R','H',...
    'A','D+N','Q+E','G','P','S'}; %One-letter abbreviation for amino acids
ratios=zeros(2,5);
xo=[sum(Var_Diet(:,PosAAPos_NHANES+39),2) sum(Var_Diet(:,NegAAPos_NHANES+39),2)];
figure;
for i=1:2
    x=xo(:,i)-[Var_Demo Var_LifeStyle ones(size(Var_Demo,1),1)]*...
        ([Var_Demo Var_LifeStyle ones(size(Var_Demo,1),1)]\xo(:,i));
    x=x(goodpos);
    q4=quantile(x,0:1/nq:1);
    label_quantile=zeros(length(goodpos),1);
    for k=1:nq
        qpos=find(x>=q4(k) & x<q4(k+1));
        if k==nq
            qpos=find(x>=q4(k));
        end
        label_quantile(qpos)=k;
        ratios(i,k)=sum(DiseaseScores(goodpos(qpos),3)==1)/length(qpos);
    end
    [~,~,p]=crosstab(label_obesity,label_quantile);
    subplot(1,2,i);
    bar(ratios(i,:),'FaceColor',[142 169 219]/255,'EdgeColor','none');
    xlabel('Dietary intake quantile');ylabel('Obesity incidence');
    ylim([0 0.6]);
    if i==1
        title(sprintf('%s\n(%s)','Total pro-obesity amino acids',...
            cell2mat(join(AACodes_NHANES(PosAAPos_NHANES),'+'))));
    else
        title(sprintf('%s\n(%s)','Total anti-obesity amino acids',...
            cell2mat(join(AACodes_NHANES(NegAAPos_NHANES),'+'))));
    end
    annotation(gcf,'textbox',[0.5*i-0.35 0.6 0.3 0.3],'String',...
        sprintf('Chi-squared\n p-value=%0.2e',p),'FitBoxToText','on','EdgeColor','none');
end

clear label_quantile qpos x xo q4 p k

%% Construct Pareto surface based on AA optimization goal derived from the NHANES data
% Anti-obesity AAs: Tryptophan, Phenylalanine, Valine,
% Aspartate+Asparagine, Glutamate+Glutamine
% Pro-obesity AAs: Lysine, Methionine, Histidine, Alanine, Glycine
[~,NegAAPos_USDA]=ismember(NegAANames,AANames);
[~,PosAAPos_USDA]=ismember(PosAANames,AANames);
c_NegAA=sum(FoodMatrix(CompletePos,AAPos(NegAAPos_USDA))/100,2);
c_PosAA=sum(FoodMatrix(CompletePos,AAPos(PosAAPos_USDA))/100,2);
DietNames={'Mediterranean','Japanese','Vegetarian','Plant-based','DASH','Paleo','Ketogenic',...
    'Atkins','American','USDA'};

NegAARange=[0 0]; %Feasible range for total intake of anti-obesity AAs
n=100; %Number of bins used to generate the feasible region
figure;

PosNegAA_NHANES=[sum(AA_Rec(:,NegAAPos_NHANES),2) sum(AA_Rec(:,PosAAPos_NHANES),2)]; %Intakes of good/bad AAs in NHANES
PosNegAA_NHANES=PosNegAA_NHANES(sum(outlier_tag,2)==0 & Age>=20,:);
Dis2Pareto=zeros(size(PosNegAA_NHANES,1),10);
colors=brewermap(10,'Set3');

list_Pareto_diets = cell(1,10);
for i_diet=1:10 % Perform the calculation for the i-th diet
    DietNames(i_diet)
    
    Diet=DietList{i_diet};
    prob_aa_diet.blx=zeros(2335,1);
    prob_aa_diet.bux=1000*ones(2335,1);
    %No more than 1kg for one food
    prob_aa_diet.blc=[Diet.ConstraintLBs;0;-Inf]; 
    prob_aa_diet.buc=[Diet.ConstraintUBs;3000;0];  
    %No more than 3kg total food, no baby food
    prob_aa_diet.a=[Diet.ConstraintMatrix(:,CompletePos);ones(1,length(CompletePos));label_baby_food];
    
    if i_diet<9
        prob_aa_diet.blc=[prob_aa_diet.blc;1800];
        prob_aa_diet.buc=[prob_aa_diet.buc;2200];
        prob_aa_diet.a=[prob_aa_diet.a;FoodMatrix(CompletePos,143)'/100];
    end
    
    prob_aa_diet.c=c_NegAA;
    [~,res]=mosekopt('maximize echo(0)',prob_aa_diet);
    NegAARange(2)=res.sol.bas.pobjval;
    [~,res]=mosekopt('minimize echo(0)',prob_aa_diet);
    NegAARange(1)=res.sol.bas.pobjval;
    xarray=NegAARange(1):(NegAARange(2)-NegAARange(1))/n:NegAARange(2);
    prob_pareto=prob_aa_diet;
    prob_pareto.a=[prob_aa_diet.a;c_NegAA'];
    prob_pareto.c=c_PosAA;
    yarray_min=zeros(size(xarray));
    yarray_max=zeros(size(xarray));
    
    Pareto_diets = zeros(2335,n+1);
    for j=1:n+1
        prob_pareto.blc=[prob_aa_diet.blc;xarray(j)-1e-6];
        prob_pareto.buc=[prob_aa_diet.buc;xarray(j)+1e-6];
        [~,res]=mosekopt('maximize echo(0)',prob_pareto);
        yarray_max(j)=res.sol.bas.pobjval;
        [~,res]=mosekopt('minimize echo(0)',prob_pareto);
        yarray_min(j)=res.sol.bas.pobjval;
        Pareto_diets(:,j) = res.sol.bas.xx;
    end
    list_Pareto_diets{i_diet} = Pareto_diets;
    subplot(2,5,i_diet);
    fill([xarray xarray(end:-1:1)],[yarray_min yarray_max(end:-1:1)],[0.9 0.9 0.9]);
    hold on;
    plot(xarray,yarray_min,'LineWidth',2,'Color',colors(6,:));
    hold on;
    scatter(PosNegAA_NHANES(:,1),PosNegAA_NHANES(:,2),[],colors(5,:),'Marker','.');%[143 170 220]/255);
    Dis2Pareto(:,i_diet)=min(pdist2([xarray(:) yarray_min(:)],PosNegAA_NHANES))';
    xlabel('AAs-to-maximize [g]');
    ylabel('AAs-to-minimize [g]');
    xlim([0 400]);ylim([0 300]);
    title(DietNames(i_diet));
end

clear xarray yarray_min yarray_max i_diet colors res prob_aa_diet prob_pareto ...
    c_NegAA c_PosAA n

% Correlate deviation from Pareto surface to obesity incidence
ratios=zeros(10,5);
p=zeros(10,1);
energy_in=Mat_NutTotal_NHANES(:,33);
energy_in=energy_in(sum(outlier_tag,2)==0 & Age>=20);
figure;
map=brewermap(64,'Blues');
for i=1:10
    f=polyfit(energy_in,Dis2Pareto(:,i),6);
    dr=Dis2Pareto(:,i)-polyval(f,energy_in);
    subplot(2,5,i);
    dscatter(energy_in,Dis2Pareto(:,i));
    xx=0:quantile(energy_in,0.999)/100:quantile(energy_in,0.999);
    hold on;
    plot(xx,polyval(f,xx));
    xlabel('Protein intake [g/day]');ylabel('Distance (D(x,PS))');title(DietNames(i));
    xlim([0 400]);ylim([0 40]);colormap(map);
    
    label_quantile=label_obesity;
    dr=dr(goodpos);
    q4=quantile(dr,0:0.2:1);
    for j=1:5
        label_quantile(dr>=q4(j) & dr<=q4(j+1))=j;
    end
    for j=1:5
        ratios(i,j)=sum(label_obesity(label_quantile==j))/sum(label_quantile==j);
    end
    [~,~,p(i)]=crosstab(label_obesity,label_quantile);
end
figure;
cmap_now = brewermap(100,'RdYlGn');
colormap(cmap_now(75:-1:26,:));
heatmap_cluster(ratios,DietNames,1:5,[0.3 0.44]);
xlabel(sprintf('Quantile of deviation\n from Pareto surface'));
title(sprintf('Association between\n obesity incidence and\n deviation from Pareto\n surface for diet'));
xtickangle(0);

%% Analyze foods appearing in different Pareto-diets
label_foods_in_Pareto_diets = zeros(2335,10);
for i=1:10
    x = list_Pareto_diets{i};
    label_foods_in_Pareto_diets(max(x,[],2)>0,i) = 1;
end

%% Save variables to the workspace file in ../data
save ../data/AA_variables.mat;
