Demo_Reduced = Demo(sum(outlier_tag,2)==0 & Age>=20,:);
sample_weights = array2table(Demo_Reduced.wtmec2yr);
independent_variables = array2table([Var_Demo Var_LifeStyle Var_Diet]);
disease_variables = array2table(DiseaseScores);

sample_weights.Properties.RowNames = string(pID_keep(sum(outlier_tag,2)==0 & Age>=20));
sample_weights.Properties.VariableNames = {'Survey weight'};

independent_variables.Properties.RowNames = string(pID_keep(sum(outlier_tag,2)==0 & Age>=20,:));
independent_variables.Properties.VariableNames = [{'Batch','Age','Gender','Race_Mexican','Race_Other_Hispanic','Race_Non_Hispanic_White',...
    'Race_Non_Hispanic_Black','Race_Other','Income','Education','Marital_Status','Smoking','Alcohol','Physical_Activity',...
    'Insurance_Coverage', 'Access_to_Healthcare'} NutAlias_NHANES AANames_NHANES];
    
disease_variables.Properties.RowNames = string(pID_keep(sum(outlier_tag,2)==0 & Age>=20,:));
disease_variables.Properties.VariableNames = DiseaseNames;

writetable(disease_variables, 'disease_variables_for_lr.csv');
writetable(independent_variables, 'independent_variables_for_lr.csv');
