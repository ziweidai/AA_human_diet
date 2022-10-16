var_importance = readtable("logistic_regression_feature_importance.csv",...
    'ReadRowNames',true);
nut_annotation = readtable("NHANES_Nutrients_Annotation.csv",'ReadRowNames',true);
nut_types = {'Energy', 'Macronutrients', 'Macronutrient subtypes',...
    'Vitamins', 'Minerals', 'Others'};

n_var = size(var_importance,1);
n_confounder = 16;
n_aa = 18;
var_type = cell(1, size(var_importance,1));
for i = 1:n_confounder
    var_type{i} = 'Demographic and life-style';
end
for i = 1:size(nut_annotation,2)
    idx = find(nut_annotation{:,i}==1) + n_confounder;
    var_type(idx) = nut_types(i);
end
for i = (n_var-n_aa+1):n_var
    var_type{i} = 'Amino acids';
end

%% Compute percentage of variables with non-zero coefficients
var_type_count = [n_confounder sum(table2array(nut_annotation)) n_aa];
var_accum_matrix = [ones(n_confounder,1) zeros(n_confounder,7);...
    zeros(n_var-n_confounder-n_aa,1) table2array(nut_annotation) zeros(n_var-n_confounder-n_aa,1);
    zeros(n_aa,7) ones(n_aa,1)];
var_type_nonzero_count = (table2array(var_importance)~=0)'*var_accum_matrix;
non_zero_percentage = var_type_nonzero_count./var_type_count;
figure;
heatmap_cluster(non_zero_percentage, var_importance.Properties.VariableNames,...
    ['Demographic and life-style variables' nut_types 'Amino acids'], [0 1]);

%% Compare variable importance between types of variables
figure;
for i = 1:4
    subplot(2,2,i);
    idx = find(var_importance{:,i} ~= 0);
    %idx = intersect(find(var_importance{:,i} ~= 0),find(sum(var_accum_matrix(:,[3 4 8]),2) ~= 0));
    data = abs(var_importance{:,i});
    violinplot(log(data(idx)), var_type(idx));
    box on;
    title(var_importance.Properties.VariableNames(i));
    ylabel('Variable importance');
    xtickangle(45);
end

