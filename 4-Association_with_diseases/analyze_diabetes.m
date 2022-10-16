% Analyze the association between diabetes and amino acid intake
%------------------------------------------------------------------------
% Variables considered:
% LBXGH - Glycohemoglobin (%)  <5.7%: Normal; 5.7%~6.4%: Prediabetes;
% >=6.5%: Diabetes
% LBXGLU - fasting plasma glucose (mg/dL) <=99: Normal; 100~125:
% Prediabetes; >=126: Diabetes
% LBXGLT - Oral glucose tolerance test (mg/dL) <=139: Normal; 140~199:
% Prediabetes; >=200: Diabetes
%------------------------------------------------------------------------
% Source: https://www.niddk.nih.gov/health-information/diabetes/overview/tests-diagnosis
DiabetesVariables=table2array(ClinVars(:,{'lbxgh','lbxglu','lbxglt'}));
DiabetesVariables=DiabetesVariables(sum(outlier_tag,2)==0 & Age>=20,:);

DiabetesScores=zeros(size(DiabetesVariables));
DiabetesScores(DiabetesVariables(:,1)>=5.7,1)=1;
DiabetesScores(DiabetesVariables(:,1)>=6.5,1)=2;
DiabetesScores(DiabetesVariables(:,2)>99,2)=1;
DiabetesScores(DiabetesVariables(:,2)>125,2)=2;
DiabetesScores(DiabetesVariables(:,3)>139,3)=1;
DiabetesScores(DiabetesVariables(:,3)>199,3)=2;
DiabetesScores(isnan(DiabetesVariables))=NaN;

Diabetes_VarNames = {'A1c','Blood glucose','GTT'};

label_diabetes = find(max(DiabetesScores,[],2) == 2);
label_prediabetes = setdiff(find(min(abs(DiabetesScores-1),[],2) == 0), label_diabetes);
label_healthy = setdiff((1:size(DiabetesScores,1))',union(label_diabetes,label_prediabetes));

%% Compute the correlation coefficients, doing confounding factor correction
% Correlate the diabetes-related variables to dietary amino acid intake in
% the entire population
figure;

[c,p]=partialcorr(DiabetesVariables,Var_Diet(:,40:57),[Var_Demo Var_LifeStyle],...
    'type','Spearman','rows','pairwise');
q_vec = mafdr(p(:), 'BHFDR', true);
q = reshape(q_vec,3,18);
c(q>0.05)=0;
cmap_now = brewermap(100,'RdBu');

subplot(3,1,1);heatmap_cluster(c,Diabetes_VarNames,AANames_NHANES,[-0.1 0.1],cmap_now(end:-1:1,:));
colorbar;
title('The entire cohort');

% Correlate the diabetes-related variables to dietary amino acid intake in
% prediabetic individuals
[c,p]=partialcorr(DiabetesVariables(label_prediabetes,:),Var_Diet(label_prediabetes,40:57),...
    [Var_Demo(label_prediabetes,:) Var_LifeStyle(label_prediabetes,:)],...
    'type','Spearman','rows','pairwise');
q_vec = mafdr(p(:), 'BHFDR', true);
q = reshape(q_vec,3,18);
c(q>0.05)=0;
cmap_now = brewermap(100,'RdBu');

subplot(3,1,2);heatmap_cluster(c,Diabetes_VarNames,AANames_NHANES,[-0.1 0.1],cmap_now(end:-1:1,:));
colorbar;
title('The prediabetic subpopulation');

% Correlate the diabetes-related variables to dietary amino acid intake in
% prediabetic individuals
[c,p]=partialcorr(DiabetesVariables(label_diabetes,:),Var_Diet(label_diabetes,40:57),...
    [Var_Demo(label_diabetes,:) Var_LifeStyle(label_diabetes,:)],...
    'type','Spearman','rows','pairwise');
q_vec = mafdr(p(:), 'BHFDR', true);
q = reshape(q_vec,3,18);
c(q>0.05)=0;
cmap_now = brewermap(100,'RdBu');

subplot(3,1,3);heatmap_cluster(c,Diabetes_VarNames,AANames_NHANES,[-0.1 0.1],cmap_now(end:-1:1,:));
colorbar;
title('The diabetic subpopulation');

%% Compute the correlation coefficients, without confounding factor correction
% Correlate the diabetes-related variables to dietary amino acid intake in
% the entire population
corr_coefs_comb = [];
cohort_names = {'The entire cohort','Prediabetic','Diabetic'};

figure;

[c,p]=corr(DiabetesVariables,Var_Diet(:,40:57),'type','Spearman','rows','pairwise');
q_vec = mafdr(p(:), 'BHFDR', true);
q = reshape(q_vec,3,18);
c(q>0.05)=0;
corr_coefs_comb = [corr_coefs_comb c'];
cmap_now = brewermap(100,'RdBu');

% Correlate the diabetes-related variables to dietary amino acid intake in
% prediabetic individuals
[c,p]=corr(DiabetesVariables(label_prediabetes,:),Var_Diet(label_prediabetes,40:57),...
    'type','Spearman','rows','pairwise');
q_vec = mafdr(p(:), 'BHFDR', true);
q = reshape(q_vec,3,18);
c(q>0.05)=0;
corr_coefs_comb = [corr_coefs_comb c'];


% Correlate the diabetes-related variables to dietary amino acid intake in
% prediabetic individuals
[c,p]=corr(DiabetesVariables(label_diabetes,:),Var_Diet(label_diabetes,40:57),...
     'type','Spearman','rows','pairwise');
q_vec = mafdr(p(:), 'BHFDR', true);
q = reshape(q_vec,3,18);
c(q>0.05)=0;
corr_coefs_comb = [corr_coefs_comb c'];

figure;
for i=1:3
    subplot(1,3,i);
    heatmap_cluster(corr_coefs_comb(:,i:3:i+6),AANames_NHANES,cohort_names,[-0.1 0.1],cmap_now(end:-1:1,:));
    title(Diabetes_VarNames(i));
end
colorbar;


%% Correlate the 2-AAs combinations with the diabetes-related variables
AA_comb_names = cell(1,9*17);
AA_combs = zeros(size(Var_Diet,1), 9*17);
AACodes_NHANES={'W','T','I','L','K','M','C','F','Y','V','R','H',...
    'A','D+N','Q+E','G','P','S'}; %One-letter abbreviation for amino acids
count = 0;
for i = 1:17
    for j = i+1:18
        count = count + 1;
        AA_comb_names(count) = join(AACodes_NHANES([i j]),'+');
        AA_combs(:,count) = Var_Diet(:,39+i) + Var_Diet(:,39+j);
    end
end

[c, p] = partialcorr(DiabetesVariables,AA_combs,[Var_Demo Var_LifeStyle],'type','Spearman','rows','pairwise');
q_vec = mafdr(p(:), 'BHFDR', true);
q = reshape(q_vec,3,9*17);     
c(q>0.05) = 0;
idx_kept = find(max(abs(c))>0);

figure;
heatmap_cluster(c(:,idx_kept),Diabetes_VarNames,AA_comb_names(idx_kept),[-0.1 0.1], cmap_now(end:-1:1,:));
colorbar;
title('The entire cohort');

[c, p] = partialcorr(DiabetesVariables(label_prediabetes,:),AA_combs(label_prediabetes,:),...
    [Var_Demo(label_prediabetes,:) Var_LifeStyle(label_prediabetes,:)],'type','Spearman','rows','pairwise');
q_vec = mafdr(p(:), 'BHFDR', true);
q = reshape(q_vec,3,9*17);     
c(q>0.05) = 0;
idx_kept = find(max(abs(c))>0);

if length(idx_kept)>1
    figure;
    heatmap_cluster(c(:,idx_kept),Diabetes_VarNames,AA_comb_names(idx_kept),[-0.1 0.1], cmap_now(end:-1:1,:));
    colorbar;
    title('The prediabetic subpopulation');
end

[c, p] = partialcorr(DiabetesVariables(label_diabetes,:),AA_combs(label_diabetes,:),...
    [Var_Demo(label_diabetes,:) Var_LifeStyle(label_diabetes,:)],'type','Spearman','rows','pairwise');
q_vec = mafdr(p(:), 'BHFDR', true);
q = reshape(q_vec,3,9*17);     
c(q>0.05) = 0;
idx_kept = find(max(abs(c))>0);

if length(idx_kept)>1
    figure;
    heatmap_cluster(c(:,idx_kept),Diabetes_VarNames,AA_comb_names(idx_kept),[-0.1 0.1], cmap_now(end:-1:1,:));
    colorbar;
    title('The diabetic subpopulation');
end
        
        
%% Correlate diabetes with other known variables that are diabetes-related
total_calories = Var_Diet(:,34);
total_carb = Var_Diet(:,9);
carb_types = Var_Diet(:,[11 21]);
bmi = BMI(sum(outlier_tag,2)==0 & Age>=20);
physical_activity = Var_LifeStyle(:,3);

control_var_names = {'Total calories','Total carbohydrate','Dietary fiber','Sugar','BMI','Physical activity'};
positive_control_variables = [total_calories total_carb carb_types bmi physical_activity];

corr_coefs_comb = [];
[c,p] = corr(positive_control_variables,DiabetesVariables,'type','Spearman','rows','pairwise');
q_vec = mafdr(p(:), 'BHFDR', true);
q = reshape(q_vec,size(c));     
c(q>0.05) = 0;
corr_coefs_comb = [corr_coefs_comb c];

[c,p] = corr(positive_control_variables(label_prediabetes,:),DiabetesVariables(label_prediabetes,:),'type','Spearman','rows','pairwise');
q_vec = mafdr(p(:), 'BHFDR', true);
q = reshape(q_vec,size(c));     
c(q>0.05) = 0;
corr_coefs_comb = [corr_coefs_comb c];

[c,p] = corr(positive_control_variables(label_diabetes,:),DiabetesVariables(label_diabetes,:),'type','Spearman','rows','pairwise');
q_vec = mafdr(p(:), 'BHFDR', true);
q = reshape(q_vec,size(c));     
c(q>0.05) = 0;
corr_coefs_comb = [corr_coefs_comb c];

figure;
for i=1:3
    subplot(1,3,i);
    heatmap_cluster(corr_coefs_comb(:,i:3:i+6),control_var_names,cohort_names,[-0.2 0.2],cmap_now(end:-1:1,:));
    title(Diabetes_VarNames(i));
end


%% Correlate the diabetes-related variables with intake of all other nutritional variables
corr_coefs_comb = [];
cohort_names = {'The entire cohort','Prediabetic','Diabetic'};

[c,p] = corr(Var_Diet(:,1:39),DiabetesVariables,'type','Spearman','rows','pairwise');
q_vec = mafdr(p(:), 'BHFDR', true);
q = reshape(q_vec,size(c));     
c(q>0.05) = 0;
corr_coefs_comb = [corr_coefs_comb c];

[c,p] = corr(Var_Diet(label_prediabetes,1:39),DiabetesVariables(label_prediabetes,:),'type','Spearman','rows','pairwise');
q_vec = mafdr(p(:), 'BHFDR', true);
q = reshape(q_vec,size(c));     
c(q>0.05) = 0;
corr_coefs_comb = [corr_coefs_comb c];

[c,p] = corr(Var_Diet(label_diabetes,1:39),DiabetesVariables(label_diabetes,:),'type','Spearman','rows','pairwise');
q_vec = mafdr(p(:), 'BHFDR', true);
q = reshape(q_vec,size(c));     
c(q>0.05) = 0;
corr_coefs_comb = [corr_coefs_comb c];

figure;
for i=1:3
    subplot(1,3,i);
    heatmap_cluster(corr_coefs_comb(:,i:3:i+6),NutAlias_NHANES,cohort_names,[-0.2 0.2],cmap_now(end:-1:1,:));
    title(Diabetes_VarNames(i));
end

%% Draw scatter plot and violin plot to compare age and the diabetes-related variables
diabetes_labels_merged = repmat({'healthy'},size(DiabetesScores,1),1);
diabetes_labels_merged(label_prediabetes) = {'prediabetic'};
diabetes_labels_merged(label_diabetes) = {'diabetic'};

figure;
for i=1:3
    subplot(1,3,i);
    gscatter(Var_Demo(:,2),DiabetesVariables(:,i),diabetes_labels_merged);
    xlabel('Age');
    ylabel(Diabetes_VarNames(i));
    title(Diabetes_VarNames(i));
end

gender_labels = repmat({'Female'},size(DiabetesScores,1),1);
gender_labels(Var_Demo(:,3)==1) = {'Male'};
figure;
cohort_labels = {1:length(gender_labels),label_prediabetes,label_diabetes};
cohort_names = {'The entire cohort','Prediabetic','Diabetic'};
for i = 1:3
    for j = 1:3
        subplot(3,3,(i-1)*3+j);
        violinplot(DiabetesVariables(cohort_labels{j},i),gender_labels(cohort_labels{j}),'ShowData',false);
        title(sprintf('%s,%s',Diabetes_VarNames{i},cohort_names{j}));
        ylabel(Diabetes_VarNames(i));
    end
end

%% Draw scatter plots comparing dietary amino acid intake to diabetes-related variables
for n = 1:3
    figure;
    for i = 1:18
        subplot(3,6,i);
        gscatter(Var_Diet(:,39+i),DiabetesVariables(:,n),diabetes_labels_merged);
        xlabel(AANames_NHANES(i));
        ylabel(Diabetes_VarNames(n));
    end
end