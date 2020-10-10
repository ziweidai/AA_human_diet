%--------------------------------------------------------------------------
% Analyze human dietary AA intake values obtained by imputation of NHANES
% data
%--------------------------------------------------------------------------
% **Notice**: Run 'Load_processed_NHANES_datasets.m' to load and process 
% the required datasets before running this script
%--------------------------------------------------------------------------


%% Validate imputation of dietary AA intakes using omics datasets
% Validation 1: comparison with daily protein intake values in NHANES
figure;
scatter(sum(AA_Rec,2),Mat_NutFood_NHANES(:,33),'Marker','.');
hold on;
plot([0 600],[0 600]);
xlabel('Daily AA intake [g]');
ylabel('Daily protein intake [g]');

% Validation 2: comparison with AA uptake fluxes of NCI-60 cancer cells and
% human plasma concentrations of AAs in HMDB
AA_NCI60=readtable('Datasets_for_validation/AA_Exchange_NCI60.csv','ReadRowNames',true);
mapAA_NCI60_to_NHANES=zeros(18,18);
for i=1:18
    mapAA_NCI60_to_NHANES(i,:)=contains(AANames_NHANES,AA_NCI60.Properties.RowNames(i));
end
mean_NCI60=nanmean(table2array(AA_NCI60)'*mapAA_NCI60_to_NHANES);
std_NCI60=nanstd(table2array(AA_NCI60)'*mapAA_NCI60_to_NHANES);

AA_HMDB=readtable('Datasets_for_validation/AA_HMDB_Weight.csv');
mapAA_HMDB_to_NHANES=zeros(20,18);
for i=1:20
    mapAA_HMDB_to_NHANES(i,:)=contains(AANames_NHANES,AA_HMDB.Properties.VariableNames(i));
end
mean_HMDB=nanmean(table2array(AA_HMDB)*mapAA_HMDB_to_NHANES);
std_HMDB=nanstd(table2array(AA_HMDB)*mapAA_HMDB_to_NHANES);

mean_NHANES=mean(AA_Prot_Ratio_Rec);
std_NHANES=std(AA_Prot_Ratio_Rec);

TitleList={'Dietary intake [g/g total AA]','Blood concentration [mg/L]','Uptake flux [fmol/cell/h]'};
MeanList={mean_NHANES,mean_HMDB,-mean_NCI60};
StdList={std_NHANES,std_HMDB,std_NCI60};
count=0;
figure;
for i=1:2
    for j=i+1:3
        xmean=MeanList{i};
        ymean=MeanList{j};
        xstd=StdList{i};
        ystd=StdList{j};
        count=count+1;
        subplot(1,3,count);
        errorbar(xmean(xmean>0 & ymean>0),ymean(xmean>0 & ymean>0),...
            ystd(xmean>0 & ymean>0),ystd(xmean>0 & ymean>0),...
            xstd(xmean>0 & ymean>0),xstd(xmean>0 & ymean>0),'o');
        xlabel(TitleList{i});
        ylabel(TitleList{j});
        [c,p]=corr(xmean(xmean>0 & ymean>0)',ymean(xmean>0 & ymean>0)',...
            'type','Spearman');
        title(sprintf('Spearman correlation=%0.2f,p=%0.2f',c,p));
    end
end

% Validation 3: compare dietary AA intake with composition of cell culture medium
AA_Medium=readtable('Datasets_for_Validation/AA_culture_medium.csv','ReadRowNames',true);
mapAA_Medium_to_NHANES=zeros(20,18);
for i=1:20
    mapAA_Medium_to_NHANES(i,:)=contains(AANames_NHANES,AA_Medium.Properties.RowNames(i));
end
MediumNames=AA_Medium.Properties.VariableNames;
AA_Medium_NHANES=mapAA_Medium_to_NHANES'*table2array(AA_Medium);
figure;
for i=1:7
    subplot(2,4,i);
    errorbar(AA_Medium_NHANES(:,i),mean_NHANES',std_NHANES','o');
    xlabel('Medium concentration [g/L]');
    ylabel('Dietary intake [g/g total AA]');
    [c,p]=corr(AA_Medium_NHANES(:,i),mean_NHANES',...
            'type','Spearman');
    title(sprintf('%s\nSpearman correlation=%0.2f,p=%0.2f',...
        replace(MediumNames{i},'_','-'),c,p));
end
clear AA_HMDB AA_Medium AA_Medium_NHANES AA_NCI60 c p count i j mapAA_HMDB_to_NHANES ...
mapAA_Medium_to_NHANES mapAA_NCI60_to_NHANES mean_HMDB mean_NCI60 mean_NHANES ...
MeanList MediumNames std_HMDB std_NCI60 std_NHANES StdList TitleList xmean xstd ...
ymean ystd

%% Plot distributions of dietary AA compositions
figure;
[~,idx]=sort(median(AA_Prot_Ratio_Rec));
violinplot(AA_Prot_Ratio_Rec(:,idx),AANames_NHANES(idx),'ShowData',false);
ylabel('Dietary intake [g/g total AA]');
xtickangle(45);xlim([0 19]);box on;
title('Distributions of dietary AA composition in NHANES');

clear idx

%% tSNE analysis of dietary AA compositions in NHANES
rng(2019); %set a fixed random number seed to generate reproducible results
x=tsne(AA_Prot_Ratio_Rec);
figure;
subplot(2,2,1);
map=brewermap(64,'RdBu');
scatter(x(:,1),x(:,2),[],Age,'Marker','.');colorbar;xlim([-100 100]);
title('Age');xlabel('tSNE 1');ylabel('tSNE 2');
colormap(map(end:-1:1,:));
set(gca,'XColor', 'none','YColor','none')
subplot(2,2,2);
map=brewermap(4,'Set3');
colors=zeros(30899,3);
for i=1:4
    colors(Batch==i,:)=repmat(map(i,:),sum(Batch==i),1);
end
rp=randperm(30899);
scatter(x(rp,1),x(rp,2),[],colors(rp,:),'Marker','.');xlim([-100 100]);
title('Batch');xlabel('tSNE 1');ylabel('tSNE 2');
set(gca,'XColor', 'none','YColor','none')
subplot(2,2,3);
map=brewermap(2,'Pastel1');
for i=1:2
    colors(Gender==i,:)=repmat(map(i,:),sum(Gender==i),1);
end
scatter(x(:,1),x(:,2),[],colors,'Marker','.');xlim([-100 100]);
title('Sex');xlabel('tSNE 1');ylabel('tSNE 2');
set(gca,'XColor', 'none','YColor','none')
subplot(2,2,4);
map=brewermap(5,'Set3');
for i=1:5
    colors(table2array(Demo(:,'ridreth1'))==i,:)=repmat(map(i,:),sum(table2array(Demo(:,'ridreth1'))==i),1);
end
scatter(x(:,1),x(:,2),[],colors,'Marker','.');xlim([-100 100]);
title('Ethnicity');xlabel('tSNE 1');ylabel('tSNE 2');
set(gca,'XColor', 'none','YColor','none')

clear i x rp colors map

%% PCA analysis of dietary AA compositions in NHANES
[~,x,~,~,explained,~]=pca(AA_Prot_Ratio_Rec);
acc_explained=zeros(1,18);
for i=1:18
    acc_explained(i)=sum(explained(1:i));
end
figure;
bar(explained,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0 0 0]);
hold on;
color=[0.5529 0.8275 0.7804];
plot(1:18,acc_explained,'-o','Color',color,'MarkerFaceColor',color);
xlabel('Number of PCs');ylabel('%Variance explained');
title('PCA of dietary AA composition in NHANES');
legend('Individual explained variance','Cumulated explained variance');

figure;
xt=strcat('PC1(',sprintf('%0.1f',explained(1)),'%)');
yt=strcat('PC2(',sprintf('%0.1f',explained(2)),'%)');
subplot(2,2,1);
map=brewermap(64,'RdBu');
scatter(x(:,1),x(:,2),[],Age,'Marker','.');colorbar;
title('Age');xlabel(xt);ylabel(yt);xticks([]);yticks([]);box on;
colormap(map(end:-1:1,:));
subplot(2,2,2);
map=brewermap(4,'Set3');
colors=zeros(30899,3);
for i=1:4
    colors(Batch==i,:)=repmat(map(i,:),sum(Batch==i),1);
end
rp=randperm(30899);
scatter(x(rp,1),x(rp,2),[],colors(rp,:),'Marker','.');
title('Batch');xlabel(xt);ylabel(yt);xticks([]);yticks([]);box on;
subplot(2,2,3);
map=brewermap(2,'Pastel1');
for i=1:2
    colors(Gender==i,:)=repmat(map(i,:),sum(Gender==i),1);
end
scatter(x(:,1),x(:,2),[],colors,'Marker','.');
title('Sex');xlabel(xt);ylabel(yt);xticks([]);yticks([]);box on;
subplot(2,2,4);
map=brewermap(5,'Set3');
for i=1:5
    colors(table2array(Demo(:,'ridreth1'))==i,:)=repmat(map(i,:),sum(table2array(Demo(:,'ridreth1'))==i),1);
end
scatter(x(:,1),x(:,2),[],colors,'Marker','.');
title('Ethnicity');xlabel(xt);ylabel(yt);xticks([]);yticks([]);box on;

clear i x xt yt map color colors explained acc_explained rp

%% Plot trends of dietary amino acid composition against age
figure;
for i=1:18
    y=AA_Prot_Ratio_Rec(:,i);    
    f=fit(Age,y,'poly3');
    subplot(3,6,i);
    scatter(Age,y,'.','MarkerEdgeColor',[0.8 0.8 0.8]);
    hold on;
    plot(0:80,f(0:80),'Color',[0 0 0],'LineWidth',2);
    title(AANames_NHANES(i));xlim([0 80]);box on;
    xlabel('Age [year]');ylabel('AA/Total AA');
end

clear i y f
    