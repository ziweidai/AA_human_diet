% Do ANOVA on different combinations of amino acids
branched_chain_AAs = {'Valine', 'Leucine', 'Isoleucine'};
essential_AAs = {'Valine', 'Leucine', 'Isoleucine', 'Phenylalanine',...
    'Methionine', 'Histidine', 'Threonine', 'Tryptophan', 'Lysine'};
ketogenic_AAs = {'Leucine', 'Lysine'};
ketogenic_glucogenic_AAs = {'Phenylalanine', 'Isoleucine', 'Threonine',...
    'Tryptophan', 'Tyrosine'};
glucogenic_AAs = setdiff(AANames,union(ketogenic_AAs,ketogenic_glucogenic_AAs));

AA_combinations = {branched_chain_AAs, essential_AAs, ketogenic_AAs,...
    ketogenic_glucogenic_AAs, glucogenic_AAs};
AA_combination_names = {'Branched-chain amino acids', 'Essential amino acids',...
    'Ketogenic amino acids', 'Amino acids that are both ketogenic and glucogenic',...
    'Glucogenic amino acids'};

F_ANOVA_AA_Comb = zeros(length(AA_combination_names),1);
for i = 1:length(AA_combination_names)
    AA_comb = sum(gAA_gProt(idx_keep,ismember(AANames,AA_combinations{i})),2);
    [~,tbl]=anova1(AA_comb,FoodTypes_Reduced(idx_keep),'off');
    F_ANOVA_AA_Comb(i) = tbl{2,5};
end
figure;
subplot(1,2,1);
bar(1:length(AA_combination_names),F_ANOVA_AA_Comb);
xticks(1:length(AA_combination_names));
xticklabels(AA_combination_names);
ylabel('F-statistic from one-way ANOVA');
title('One-way ANOVA comparing AA subtypes across foods');

%% Perform additional ANOVA analysis on the 500,000 randomly sampled diets
F_ANOVA_AA_Comb_Diets = zeros(length(AA_combination_names),1);
aa_ratio_samp_comb = [];
for i = 1:10
    aa_ratio_samp = list_aa_samp{i};
    aa_ratio_samp = aa_ratio_samp./sum(aa_ratio_samp,2);
    aa_ratio_samp_comb = [aa_ratio_samp_comb;aa_ratio_samp];
end
for i = 1:length(AA_combination_names)
    AA_comb = sum(aa_ratio_samp_comb(:,ismember(AANames,AA_combinations{i})),2);
    [~,tbl] = anova1(reshape(AA_comb,50000,10),1:10,'off');
    F_ANOVA_AA_Comb_Diets(i) = tbl{2,5};
end
subplot(1,2,2);
bar(1:length(AA_combination_names),F_ANOVA_AA_Comb_Diets);
xticks(1:length(AA_combination_names));
xticklabels(AA_combination_names);
ylabel('F-statistic from one-way ANOVA');
title('One-way ANOVA comparing AA subtypes across diets');

%% Perform ANOVA analysis on randomly selected three dietary patterns
F_ANOVA_AA_Random_3_Diets = zeros(10,18);
for i = 1:10
    random_diets = randperm(10,3);
    for j = 1:18
        aa_comb = reshape(aa_ratio_samp_comb(:,j),50000,10);
        [~,tbl] = anova1(aa_comb(:,random_diets),1:3,'off');
        F_ANOVA_AA_Random_3_Diets(i,j) = tbl{2,5};
    end
end

%% Analyze overlap between vitamin D-rich and protein-rich foods
protein_in_foods = FoodMatrix(CompletePos,140);
vitamin_d_in_foods = FoodMatrix(CompletePos,171);
protein_in_NHANES = Mat_NutTotal_NHANES(:,33);
vitamin_d_in_NHANES = Mat_NutTotal_NHANES(:,38);
figure;
subplot(1,2,1);
scatter(protein_in_foods, vitamin_d_in_foods, 'Marker',...
    '.','MarkerEdgeColor','#70AD47');
title('Correlation of protein and vitamin D levels in foods');
xlabel('Protein (g/100g food)');ylabel('Vitamin D (\mu g)/100g food');
box on;
subplot(1,2,2);
scatter(protein_in_NHANES,vitamin_d_in_NHANES, 'Marker',...
    '.','MarkerEdgeColor','#70AD47');
title('Correlation of protein and vitamin D intake in individuals');
xlabel('Protein (g/100g food)');ylabel('Vitamin D (\mu g)/100g food');
box on;

%% Get the loadings from PCA of AA profiles in foods
[coeff,score,~,~,explained,~]=pca(gAA_gProt(Foods_wProt,:));

%% Perform ANOVA analysis with the log-transformed carb, fat, and AA values in foods
F_ANOVA_AA_log=zeros(18,1);
F_ANOVA_Carb_log=zeros(3,1);
F_ANOVA_Fat_log=zeros(9,1);
for i=1:18
    [~,tbl]=anova1(log(1+gAA_gProt(idx_keep,i)),FoodTypes_Reduced(idx_keep),'off');
    F_ANOVA_AA_log(i)=tbl{2,5};
end
for i=1:3
    [~,tbl]=anova1(log(1+CarbRatio(CompletePos(idx_keep),i)),FoodTypes_Reduced(idx_keep),'off');
    F_ANOVA_Carb_log(i)=tbl{2,5};
end
for i=1:9
    [~,tbl]=anova1(log(1+LipidRatio_Combine(CompletePos(idx_keep),i)),FoodTypes_Reduced(idx_keep),'off');
    F_ANOVA_Fat_log(i)=tbl{2,5};
end
F_ANOVA_AA_Comb_log = zeros(length(AA_combination_names),1);
for i = 1:length(AA_combination_names)
    AA_comb = log(1+sum(gAA_gProt(idx_keep,ismember(AANames,AA_combinations{i})),2));
    [~,tbl]=anova1(AA_comb,FoodTypes_Reduced(idx_keep),'off');
    F_ANOVA_AA_Comb_log(i) = tbl{2,5};
end

%% Perform ANOVA analysis with the log-transformed carb, fat, and AA values in diets
F_ANOVA_AA_Comb_Diets_log = zeros(length(AA_combination_names),1);
F_ANOVA_AA_Diets_log = zeros(18,1);
F_ANOVA_Carb_Diet_log = zeros(3,1);
F_ANOVA_Fat_Diet_log = zeros(3,1);
aa_ratio_samp_comb = [];
for i = 1:10
    aa_ratio_samp = list_aa_samp{i};
    aa_ratio_samp = aa_ratio_samp./sum(aa_ratio_samp,2);
    aa_ratio_samp_comb = [aa_ratio_samp_comb;aa_ratio_samp];
end
for i = 1:length(AA_combination_names)
    AA_comb = log(1+sum(aa_ratio_samp_comb(:,ismember(AANames,AA_combinations{i})),2));
    [~,tbl] = anova1(reshape(AA_comb,50000,10),1:10,'off');
    F_ANOVA_AA_Comb_Diets_log(i) = tbl{2,5};
end
for i = 1:18
    [~,tbl] = anova1(reshape(log(1+aa_ratio_samp_comb(:,i)),50000,10),1:10,'off');
    F_ANOVA_AA_Diets_log(i) = tbl{2,5};
end
for i = 1:3
    [~,tbl] = anova1(log(1+list_carb_ratio{i}),1:10,'off');
    F_ANOVA_Carb_Diet_log(i) = tbl{2,5};
    [~,tbl] = anova1(log(1+list_lipid_ratio{i}),1:10,'off');
    F_ANOVA_Fat_Diet_log(i) = tbl{2,5};
end

%% Perform PCA on the log-transformed AA profiles in foods
[~,score,~,~,explained,~]=pca(log(1+gAA_gProt(Foods_wProt,:)));
figure;gscatter(score(:,1),score(:,2),CATs,map,[],20);
title('PCA of log-transformed amino acid profiles in foods');

%% PCA of log-transformed amino acid compositions of diets
diet_tag=reshape(repmat(DietNames,5000,1),50000,1);
aa_ratio_downsample=zeros(50000,18);
for i=1:10
    rp=randperm(50000);
    a=list_aa_samp{i};
    aa_ratio_downsample((i-1)*5000+1:i*5000,:)=a(rp(1:5000),:)./sum(a(rp(1:5000),:),2);    
end

[~,score,~,~,explained,~]=pca(log(1+aa_ratio_downsample));
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
title('PCA of log-transformed amino acid compositions of dietary patterns');

%% Analyze amino acid intake profiles in individuals eating ketogenic or Atkins diets
nhanes_carb_intake = Mat_NutTotal_NHANES(:,9);
nhanes_fat_intake = Mat_NutTotal_NHANES(:,13);
nhanes_energy_intake = Mat_NutTotal_NHANES(:,34);
nhanes_carb_contrb_energy = 4*nhanes_carb_intake./nhanes_energy_intake;
nhanes_fat_contrb_energy = 9*nhanes_fat_intake./nhanes_energy_intake;

dev_from_ketogenic = [nhanes_carb_contrb_energy nhanes_fat_contrb_energy] - [0.05 0.7];
dev_from_ketogenic(nhanes_carb_contrb_energy<0.05,1) = 0;
dist_nhanes_2_ketogenic_diet = (sum(dev_from_ketogenic.^2, 2)).^0.5;
corr_keto_AA_nhanes = corr(AA_Rec(Age>18,:)./sum(AA_Rec(Age>18,:),2),...
    -dist_nhanes_2_ketogenic_diet(Age>18),'type','Spearman');

aa_adult = AA_Rec(Age>18,:)./sum(AA_Rec(Age>18,:),2);
[~,idx] = sort(dist_nhanes_2_ketogenic_diet(Age>18));
aa_adult_sorted_by_ketoscore = aa_adult(idx,:);
diff_aa_keto_vs_nonketo = mean(aa_adult_sorted_by_ketoscore(1:1000,:))...
    -mean(aa_adult_sorted_by_ketoscore(1001:end,:));
diff_aa_keto_vs_nonketo = diff_aa_keto_vs_nonketo./mean(aa_adult_sorted_by_ketoscore(1:1000,:));

%% Compare AA profiles of animal- and plant-derived protein
plant_cats = {'Cereals', 'Legumes', 'Nuts', 'Vegetables'};
animal_cats = {'Beef', 'Dairy', 'Eggs', 'Lamb, Veal, and Game', 'Pork',...
    'Poultry', 'Seafood'};
plant_idx = find(ismember(FoodTypes_Reduced, plant_cats));
animal_idx = find(ismember(FoodTypes_Reduced, animal_cats));
figure;
cv_mat = zeros(3,18);
for i = 1:18
    subplot(3,6,i);
    data = [gAA_gProt(plant_idx,i);gAA_gProt(animal_idx,i)];
    violinplot(data,[repmat({'Plant'},length(plant_idx),1);...
        repmat({'Animal'},length(animal_idx),1)]);
    title(AANames(i));
    ylabel('g/g total AAs');
    box on;
    cv_mat(1,i) = std(gAA_gProt(plant_idx,i))/mean(gAA_gProt(plant_idx,i));
    cv_mat(2,i) = std(gAA_gProt(animal_idx,i))/mean(gAA_gProt(animal_idx,i));
    cv_mat(3,i) = std(data)/mean(data);
end