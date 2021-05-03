%--------------------------------------------------------------------------
% Estimate correlation between health-related variables and dietary
% composition of amino acids and other dietary variables
%--------------------------------------------------------------------------

%% Load existing variables from workspace saved in ../data
load ../data/AA_variables.mat;

%% Extract health-related variables from the NHANES dataset
% 1 - Diabetes
% 2 - Hypertension
% 3 - Obesity
% 4 - Cancer

% Diabetes
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
DiabetesScores=zeros(size(DiabetesVariables));
DiabetesScores(DiabetesVariables(:,1)>=5.7,1)=1;
DiabetesScores(DiabetesVariables(:,1)>=6.5,1)=2;
DiabetesScores(DiabetesVariables(:,2)>99,2)=1;
DiabetesScores(DiabetesVariables(:,2)>125,2)=2;
DiabetesScores(DiabetesVariables(:,3)>139,3)=1;
DiabetesScores(DiabetesVariables(:,3)>199,3)=2;
DiabetesScores(isnan(DiabetesVariables))=NaN;
% Use Glycohemoglobin only because the other variables have too many
% missing values
Diabetes=double(DiabetesScores(:,1)==2);
Diabetes(isnan(DiabetesScores(:,1)))=NaN;

clear DiabetesVariables DiabetesScores

% Hypertension
%-------------------------------------------------------------------------
% Variables considered:
% BPXSY1/2/3/4: Systolic blood pressure, 4 readings. cutoff=120mmHg
% BPXDI1/2/3/4: Diastolic blood pressure, 4 readings. cutoff=80mmHg
%-------------------------------------------------------------------------
% Source: http://www.heart.org/HEARTORG/Conditions/HighBloodPressure/GettheFactsAbou
% tHighBloodPressure/How-High-Blood-Pressure-is-Diagnosed_UCM_301873_Article.jsp#.W0zil9JKg2w
%-------------------------------------------------------------------------
BloodPressures=table2array(ClinVars(:,{'bpxsy1','bpxsy2','bpxsy3',...
    'bpxdi1','bpxdi2','bpxdi3'}));
BPTags=zeros(size(BloodPressures));
BPTags(:,1:3)=(BloodPressures(:,1:3)>120);
BPTags(:,4:6)=(BloodPressures(:,4:6)>80);
BPTags(BloodPressures*0~=0)=NaN;
Hypertension=double(sum(BPTags,2)==6);
Hypertension(isnan(sum(BPTags,2)))=NaN;
clear BPTags BloodPressures

% Obesity
%-------------------------------------------------------------------------
% Variables considered:
% BMXBMI - BMI. <18.5: underweight; 18.5~24.9: normal weight; 25.0~29.9:
% overweight; 30.0~39.9: obese; >=40: extremely obese
% Age - age of the participant.
% Additional note: for children and teens < 19 years old, the criterion for
% diagnosing obesity is defined based on the quantiles of BMI in the same
% gender-age group. <5%: underweight: 5%~85%: normal weight; 85%~95%:
% overweight; >95%: obese
%-------------------------------------------------------------------------
% Source: https://www.nichd.nih.gov/health/topics/obesity/conditioninfo/diagnosed#f1
%-------------------------------------------------------------------------
BMI=table2array(ClinVars(:,'bmxbmi'));
Obesity=zeros(size(BMI));
%Calculate obesity score for adults
Obesity(Age>=19 & BMI<18.5)=1;
Obesity(Age>=19 & BMI>=18.5 & BMI<25)=2;
Obesity(Age>=19 & BMI>=25 & BMI<30)=3;
Obesity(Age>=19 & BMI>=30 & BMI<40)=4;
Obesity(Age>=19 & BMI>=40)=5;
%Calculate obesity score for children and teens
child_groups=[repmat(0:18,1,2);ones(1,19) 2*ones(1,19)];
child_bmi_cutoffs=zeros(3,38);
for i=1:38
    cbmi_group=BMI(fix(Age)==child_groups(1,i) & Gender==child_groups(2,i));
    child_bmi_cutoffs(:,i)=quantile(cbmi_group,[0.05 0.85 0.95])';
end
child_bmi=BMI(Age<19);
child_age=Age(Age<19);
child_gender=Gender(Age<19);
child_score=zeros(size(child_bmi));
for i=1:length(child_bmi)
    cutoffs=child_bmi_cutoffs(:,fix(child_age(i))==child_groups(1,:) ...
        & child_gender(i)==child_groups(2,:));
    if child_bmi(i)<cutoffs(1)
        child_score(i)=1;
    elseif child_bmi(i)<cutoffs(2)
        child_score(i)=2;
    elseif child_bmi(i)<cutoffs(3)
        child_score(i)=3;
    else
        child_score(i)=4;
    end
end
Obesity(Age<19)=child_score;
Obesity(isnan(BMI))=NaN;
Obesity(Obesity<4)=0;
Obesity(Obesity>=4)=1;
clear cbmi_group child_age child_bmi child_bmi_cutoffs child_gender ...
    child_groups child_score cutoffs i

%Cancer (variable: MCQ220)
Cancer=table2array(ClinVars(:,'mcq220'));
Cancer(Cancer==2)=0;
Cancer(Cancer~=0 & Cancer~=1)=NaN;

%% Correct for potential confounders including demographic and life-style-related variables
% Demographic variables
Var_Demo=[Batch Age Gender Race Income Education MaritalStatus];
for i=1:size(Var_Demo,2)
    x=Var_Demo(:,i);
    xi=x;
    xi(isnan(x))=mean(x(~isnan(x)));
    Var_Demo(:,i)=xi;
end

% Life style variables: smoking status, alcohol consumption, physical activity
Var_LifeStyle=[Smoking Alcohol PhysicalActivity];
% Impute missing values using mean value
for i=1:size(Var_LifeStyle,2)
    x=Var_LifeStyle(:,i);
    xi=x;
    xi(isnan(x))=mean(x(~isnan(x)));
    Var_LifeStyle(:,i)=xi;
end
clear x xi

%Dietary variables: nutrient intakes and amino acid compositions
Var_Diet=[Mat_NutTotal_NHANES AA_Prot_Ratio_Rec];

%Detect outliers in the dietary variables
outlier_tag=zeros(size(Var_Diet));
for i=1:size(Var_Diet,2)
    x=Var_Diet(:,i);
    outlier_tag(:,i)=(x>quantile(x,0.99)*3);
end

% Remove outliers and records for non-adults (i.e. <20 years old)
Var_Diet=Var_Diet(sum(outlier_tag,2)==0 & Age>=20,:);
Var_Demo=Var_Demo(sum(outlier_tag,2)==0 & Age>=20,:);
Var_LifeStyle=Var_LifeStyle(sum(outlier_tag,2)==0 & Age>=20,:);

%Dietary variables corrected for variance explained by demographic and life
%style variables
n_left=size(Var_Diet,1);
Var_Diet_Adjusted=Var_Diet-[Var_Demo Var_LifeStyle ones(n_left,1)]*...
([Var_Demo Var_LifeStyle ones(n_left,1)]\Var_Diet);

%Create list of respond variables (cancer, diabetes, obesity, hypertension)
DiseaseNames={'Cancer','Diabetes','Obesity','Hypertension'};
DiseaseScores=[Cancer Diabetes Obesity Hypertension];
DiseaseScores=DiseaseScores(sum(outlier_tag,2)==0 & Age>=20,:);


%% Test the ability of different groups of nutrients to predict disease
% Define groups of nutritional variables
Nut_Group_Names={'Energy','Macronutrients','Macronutrient subtypes','Vitamins',...
    'Minerals','Others','Amino acid compositions'};
Nut_Groups=cell(1,7);
Nut_Groups_Adjusted=cell(1,7);
NHANES_Nut_Annotation=table2array(readtable('input_data/NHANES_Nutrients_Annotation.csv',...
    'ReadRowNames',true));
for i=1:6
    Nut_Groups{i}=Var_Diet(:,NHANES_Nut_Annotation(:,i)==1);
    Nut_Groups_Adjusted{i}=Var_Diet_Adjusted(:,NHANES_Nut_Annotation(:,i)==1);
end
Nut_Groups{7}=Var_Diet(:,40:57);
Nut_Groups_Adjusted{7}=Var_Diet_Adjusted(:,40:57); 
x=Nut_Groups{3};
y=Nut_Groups{2};
x(:,[5 6])=x(:,[5 6])./y(:,1);
x(:,[1:4 7:10])=x(:,[1:4 7:10])./y(:,2);
z=x-[Var_Demo Var_LifeStyle ones(n_left,1)]*...     %Adjust for confounders
    ([Var_Demo Var_LifeStyle ones(n_left,1)]\x);
Nut_Groups_Adjusted{2}=z;
Carb_Fat_Frac=x; %Fraction of carbohydrate and fat subtypes
Carb_Fat_Names={'PFA 18:2','Total MFA','PFA 18:3','SFA 16:0','Dietary fiber',...
    'Total sugars','Total PFA','MFA 18:1','SFA 18:0','Total SFA'};
clear x y z

% Calculate AUCs for predicting health-related variable using logistic
% regression model with dietary variables as inputs
AUCs_All=zeros(10,28);
for ni=1:10
    Nut_Group_AUCs=zeros(4,7);
    for i=1:4
        y=DiseaseScores(:,i);
        AvalPos=find(~isnan(y));
        y=y(AvalPos);
        CV_Tag=randi(5,length(AvalPos),1);
        for j=1:7
            X=Nut_Groups_Adjusted{j};
            X=X(AvalPos,:);
            ypred=zeros(length(AvalPos),1);
            for k=1:5
                test_pos=find(CV_Tag==k);
                train_pos=find(CV_Tag~=k);
                
                %{
                % Use Lasso generalized linear model (with regularization)
                [B,FitInfo]=lassoglm(X(train_pos,:),y(train_pos),'binomial','CV',3);%,'constant','off');
                idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
                B0 = FitInfo.Intercept(idxLambdaMinDeviance);
                b = [B0; B(:,idxLambdaMinDeviance)];
                %}
                
                % Use standard generalized linear model without
                % regularization
                b = glmfit(X(train_pos,:),y(train_pos),'binomial');
                
                ypred(test_pos)=glmval(b,X(test_pos,:),'logit');%,'constant','off');                
            end
            [~,~,~,Nut_Group_AUCs(i,j)]=perfcurve(y,ypred,1);            
        end
    end
    AUCs_All(ni,:)=Nut_Group_AUCs(:)';
end
mean_AUCs=reshape(mean(AUCs_All),4,7);
std_AUCs=reshape(std(AUCs_All),4,7);
clear AUCs_All y ypred Nut_Group_AUCs ni i j k n_left b X CV_Tag AvalPos ...
    train_pos test_pos

%% Compute correlation between dietary variables and health
[c,p]=partialcorr(DiseaseScores,Var_Diet(:,40:57),[Var_Demo Var_LifeStyle],...
    'type','Spearman','rows','pairwise');
c(p>0.05)=NaN;
cmap_now = brewermap(100,'RdYlGn');

figure;heatmap_cluster(c,DiseaseNames,AANames_NHANES,[-0.1 0.1],cmap_now(75:-1:26,:));
colorbar;
title('Correlation between dietary AA composition and human diseases');

[c,p]=partialcorr(DiseaseScores,Carb_Fat_Frac,[Var_Demo Var_LifeStyle],...
    'type','Spearman','rows','pairwise');
c(p>0.05)=NaN;
figure;heatmap_cluster(c,DiseaseNames,Carb_Fat_Names,[-0.1 0.1],cmap_now(75:-1:26,:));
colorbar;
title('Correlation between dietary carbohydrate and fat composition and human diseases');

clear cmap_now c

%% Save variables to the workspace file in ../data
save ../data/AA_variables.mat;