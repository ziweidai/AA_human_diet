%Construct constraints on each type of diet. The constraints are either on
%total amount of foods from one food type or on daily intake of nutrients.

%Extract serving information and calculate reciprocal of minimal servings
ReciprocalServing=1./MinServing;
ReciprocalServing(MinServing==-1)=0; %Foods without serving information were removed

% Create numeric IDs for each food type
[~,FoodTypeIDs]=ismember(FoodTypes,unique(FoodTypes));

%-----------------------------------------------------------------------
%Definition of diets
%Constraints on diets are stored in data structures named by the diet
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
%Diet 1: Mediterranean diet
%Recommended composition of a Mediterranean diet is taken from the
%Mediterranean diet pyramid today. It included recommended consumption of
%foods in several food groups.
%Source: https://www.cambridge.org/core/journals/public-health-nutrition/
%article/mediterranean-diet-pyramid-today-science-and-cultural-updates/
%70359644D12A038AC003B935AA04E120
%-----------------------------------------------------------------------
ConstraintNames={'Cereals','Fish','Fruits','Legumes','White meat',...
    'Processed meat','Sweets','Vegetables','Potatoes','Red meat',...
    'Olives/nuts/seeds','Olive oil','Eggs','Dairy','Others'};
%Lower and upper bounds of servings of foods in each category
ConstraintLBs=[3;0.29;3;0.29;0.29;0;0;6;0;0;1;1;0.29;2;0];
ConstraintUBs=[6;Inf;6;Inf;0.29;0.14;0.29;Inf;0.43;0.29;2;Inf;0.57;2;0];
%Coefficient matrix for the constraints
ConstraintMatrix=zeros(15,8788);
%Cereals
ConstraintMatrix(1,:)=ReciprocalServing'.*(FoodTypeIDs'==7);
%Fish
ConstraintMatrix(2,:)=ReciprocalServing'.*(FoodTypeIDs'==21);
%Fruits
ConstraintMatrix(3,:)=ReciprocalServing'.*(FoodTypeIDs'==12 & contains(FoodNames','juice')==0);
%Legumes
ConstraintMatrix(4,:)=ReciprocalServing'.*(FoodTypeIDs'==14);
%White meat
ConstraintMatrix(5,:)=ReciprocalServing'.*(FoodTypeIDs'==18);
%Processed meat
ConstraintMatrix(6,:)=ReciprocalServing'.*(FoodTypeIDs'==20);
%Sweets
ConstraintMatrix(7,:)=ReciprocalServing'.*(FoodTypeIDs'==25);
%Vegetables
ConstraintMatrix(8,:)=ReciprocalServing'.*(FoodTypeIDs'==26);
%Potatoes
ConstraintMatrix(9,:)=ReciprocalServing'.*(contains(FoodNames','Potatoes'));
%Red meat
ConstraintMatrix(10,:)=ReciprocalServing'.*(FoodTypeIDs'==4 | FoodTypeIDs'==13 | FoodTypeIDs'==17);
%Olives/nuts/seeds
ConstraintMatrix(11,:)=ReciprocalServing'.*(FoodTypeIDs'==16 | contains(FoodNames','Olives'));
%Olive oil
ConstraintMatrix(12,:)=ReciprocalServing'.*(contains(FoodNames','Oil, olive'));
%Eggs
ConstraintMatrix(13,:)=ReciprocalServing'.*(FoodTypeIDs'==9);
%Dairy
ConstraintMatrix(14,:)=ReciprocalServing'.*(FoodTypeIDs'==8);
%Others
ConstraintMatrix(15,:)=(max(ConstraintMatrix(1:14,:))==0);
%-------------------------------------------------------------------------
%Construct the data structure
%-------------------------------------------------------------------------
Mediterranean.ConstraintNames=ConstraintNames;
Mediterranean.ConstraintLBs=ConstraintLBs;
Mediterranean.ConstraintUBs=ConstraintUBs;
Mediterranean.ConstraintMatrix=ConstraintMatrix;
%-------------------------------------------------------------------------
%Finished definition of Mediterranean diet
%-------------------------------------------------------------------------

%Diet 2: Japanese diet
%Recommended composition of balanced Japanese diet was adopted from the
%Japanese dietary guidelines
%Source: https://www.mobiefit.com/blog/follow-a-japanese-spinning-top-diet-to-live-longer/
ConstraintNames={'Cereals','Fruits','Vegetables','Dairy','Proteins','Others'};
%Here 'Proteins' include fish, meat, egg and soy bean
ConstraintLBs=[5;2;5;2;3;0];
ConstraintUBs=[7;2;6;2;5;0];
ConstraintMatrix=zeros(6,8788);
%Cereals
ConstraintMatrix(1,:)=ReciprocalServing'.*(FoodTypeIDs'==7);
%Fruits
ConstraintMatrix(2,:)=ReciprocalServing'.*(FoodTypeIDs'==12 & contains(FoodNames','juice')==0);
%Vegetables
ConstraintMatrix(3,:)=ReciprocalServing'.*(FoodTypeIDs'==26);
%Dairy
ConstraintMatrix(4,:)=ReciprocalServing'.*(FoodTypeIDs'==8);
%Proteins
ConstraintMatrix(5,:)=ReciprocalServing'.*(ismember(FoodTypeIDs',[4 13 17 18 21]) ...
| contains(FoodNames','Soybean'));
%Others
ConstraintMatrix(6,:)=(max(ConstraintMatrix(1:5,:))==0);
%-------------------------------------------------------------------------
%Construct the data structure
%-------------------------------------------------------------------------
Japanese.ConstraintNames=ConstraintNames;
Japanese.ConstraintLBs=ConstraintLBs;
Japanese.ConstraintUBs=ConstraintUBs;
Japanese.ConstraintMatrix=ConstraintMatrix;
%-------------------------------------------------------------------------
%Finished definition of Japanese diet
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
%Diet 3: Vegetarian diet
%Vegetarian diet is defined by zero consumption of animal products
%-------------------------------------------------------------------------
ConstraintNames={'Animal products'};
ConstraintLBs=0;
ConstraintUBs=0;
ConstraintMatrix=double(ismember(FoodTypeIDs',[1 2 4 10 13 15 17 18 19 20 21 22]));
%-------------------------------------------------------------------------
%Construct the data structure
%-------------------------------------------------------------------------
Vegetarian.ConstraintNames=ConstraintNames;
Vegetarian.ConstraintLBs=ConstraintLBs;
Vegetarian.ConstraintUBs=ConstraintUBs;
Vegetarian.ConstraintMatrix=ConstraintMatrix;
%-------------------------------------------------------------------------
%Finished definition of vegetarian diet
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
%Diet 4: Paleo diet
%Paleo diet is defined by only eating what ancient human ate. These include
%game meat, fish/seafood, fresh fruit and vegetables, nuts/seeds
%and honey. Potatoes, legumes and fruit juices are avoided
%-------------------------------------------------------------------------
ConstraintNames={'Modern products'};
ConstraintLBs=0;
ConstraintUBs=0;
ConstraintMatrix=double(1-(contains(FoodNames','Game meat') ...
    | ismember(FoodTypeIDs',[12 16 21 26]) | contains(FoodNames','19296, Honey'))) ...
    | contains(FoodNames','juice') | contains(FoodNames','Potatoes');
%-------------------------------------------------------------------------
%Construct the data structure
%-------------------------------------------------------------------------
Paleo.ConstraintNames=ConstraintNames;
Paleo.ConstraintLBs=ConstraintLBs;
Paleo.ConstraintUBs=ConstraintUBs;
Paleo.ConstraintMatrix=ConstraintMatrix;
%-------------------------------------------------------------------------
%Finished definition of Paleo diet
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
%Diet 5: Ketogenic diet
%-------------------------------------------------------------------------
%In a ketogenic diet fats contribute to >70% of the total calories while
%carbohydrates contribute to <5% of the total calories
ConstraintNames={'Fat %Calories','Carbohydrate %Calories'};
ConstraintLBs=[0;-Inf];
ConstraintUBs=[Inf;0];
ConstraintMatrix=[9*FoodMatrix(:,43)'-0.7*FoodMatrix(:,143)';...
    4*FoodMatrix(:,34)'-0.05*FoodMatrix(:,143)'];
%-------------------------------------------------------------------------
%Construct the data structure
%-------------------------------------------------------------------------
Ketogenic.ConstraintNames=ConstraintNames;
Ketogenic.ConstraintLBs=ConstraintLBs;
Ketogenic.ConstraintUBs=ConstraintUBs;
Ketogenic.ConstraintMatrix=ConstraintMatrix;
%-------------------------------------------------------------------------
%Finished definition of ketogenic diet
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
%Diet 6: Atkins diet
%-------------------------------------------------------------------------
%Atkins diet allows less than 20g carbohydrates per day
%-------------------------------------------------------------------------
ConstraintNames={'Carbohydrate'};
ConstraintLBs=0;
ConstraintUBs=20;
ConstraintMatrix=0.01*FoodMatrix(:,34)';
%-------------------------------------------------------------------------
%Construct the data structure
%-------------------------------------------------------------------------
Atkins.ConstraintNames=ConstraintNames;
Atkins.ConstraintLBs=ConstraintLBs;
Atkins.ConstraintUBs=ConstraintUBs;
Atkins.ConstraintMatrix=ConstraintMatrix;
%-------------------------------------------------------------------------
%Finished definition of Atkins diet
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
%Diet 7: American diet
%-------------------------------------------------------------------------
%American diet is defined based on the National Geographic website
%Source: https://www.nationalgeographic.com/what-the-world-eats/
%-------------------------------------------------------------------------
%Calories: 3641 per day
%Meats: 381g per day
%Sugar&sweeteners: 166g per day
%Oils: 114g per day
%Grains: 291g per day
%Starchy roots: 164g per day (potatoes,sweet potatoes,yams,taro,cassava)
%(need to be excluded from the vegetable group)
%Vegetables: 310g per day
%Fruits: 266g per day
%Eggs: 38g per day
%Milk: 704g per day
%Legumes: 9g per day
%Alcoholic beverages: 258g
%Others: 28g
%------------------------------------------------------------------------
%Constraints are set based on the values above (+/-20%)
%------------------------------------------------------------------------
ConstraintNames={'Calories','Meats','Sweets','Oils','Grains',...
    'Starchy roots','Vegetables','Fruits','Eggs','Milk','Legumes',...
    'Alcoholic beverages','Others'};
Center=[3641;381;166;114;291;164;310;266;38;704;9;258;28];
ConstraintLBs=Center*0.8;
ConstraintUBs=Center*1.2;
clear Center;
ConstraintMatrix=zeros(13,8788);
%Calories
ConstraintMatrix(1,:)=FoodMatrix(:,143)'/100;
%Meats
ConstraintMatrix(2,:)=ismember(FoodTypeIDs',[4 13 17 18 20 21]);
%Sweets
ConstraintMatrix(3,:)=(FoodTypeIDs'==25);
%Oils
ConstraintMatrix(4,:)=(FoodTypeIDs'==11);
%Grains
ConstraintMatrix(5,:)=(FoodTypeIDs'==7);
%Starchy roots
ConstraintMatrix(6,:)=(contains(FoodNames','Potato') | ...
    contains(FoodNames','Sweet potato') | contains(FoodNames','Yam') | ...
    contains(FoodNames','yam') | contains(FoodNames','Taro') | ...
    contains(FoodNames','Cassava'));
%Vegetables
ConstraintMatrix(7,:)=(FoodTypeIDs'==26)-ConstraintMatrix(6,:);
ConstraintMatrix(ConstraintMatrix<0)=0;
%Fruits
ConstraintMatrix(8,:)=(FoodTypeIDs'==12);
%Eggs
ConstraintMatrix(9,:)=(FoodTypeIDs'==9);
%Milk
ConstraintMatrix(10,:)=(FoodTypeIDs'==8);
%Legumes
ConstraintMatrix(11,:)=(FoodTypeIDs'==14);
%Alcoholic beverages
ConstraintMatrix(12,:)=contains(FoodNames','Alcoholic');
%Others
ConstraintMatrix(13,:)=1-max(ConstraintMatrix(2:12,:));
%-------------------------------------------------------------------------
%Construct the data structure
%-------------------------------------------------------------------------
American.ConstraintNames=ConstraintNames;
American.ConstraintLBs=ConstraintLBs;
American.ConstraintUBs=ConstraintUBs;
American.ConstraintMatrix=ConstraintMatrix;
%-------------------------------------------------------------------------
%Finished definition of American diet
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
%Diet 8: DASH diet
%-------------------------------------------------------------------------
%The DASH diet emphasizes reduced intake of saturated/trans fats.
%Constraints that defining the DASH diet are constructed based on the DASH
%diet pyramid and the US News website. Values from the US News are directly
%used since their units are consistent with the USDA SR database
%Source1: https://www.languageofdesires.com/dash-diet/
%Source2: https://health.usnews.com/best-diet/dash-diet
%-------------------------------------------------------------------------
ConstraintNames={'Grains','Vegetables','Fruits','Low-fat dairy','Lean meat',...
    'Nuts,seeds and legumes','Fats and oils','Sweets','Sodium','Others'};
ConstraintLBs=[6;4;4;2;0;4/7;2;0;0;0];
ConstraintUBs=[8;5;5;3;6;5/7;3;5/7;2300;0];
ConstraintMatrix=zeros(10,8788);
%Grains
ConstraintMatrix(1,:)=ReciprocalServing'.*(FoodTypeIDs'==7);
%Vegetables
ConstraintMatrix(2,:)=ReciprocalServing'.*(FoodTypeIDs'==26);
%Fruits
ConstraintMatrix(3,:)=ReciprocalServing'.*(FoodTypeIDs'==12 & contains(FoodNames','juice')==0);
%Low-fat diary
ConstraintMatrix(4,:)=ReciprocalServing'.*(FoodTypeIDs'==8 ...
    & contains(FoodNames','fat') & ~contains(FoodNames','3.'));
%Lean meat (defined by protein_fat_ratio>4)
protein_fat_ratio=FoodMatrix(:,140)./FoodMatrix(:,43); %mass ratio of protein to fat
ConstraintMatrix(5,:)=ReciprocalServing'.*...
    (ismember(FoodTypeIDs',[4 13 17 18 20 21]) & protein_fat_ratio'>4);
%Nuts,seeds and legumes
ConstraintMatrix(6,:)=ReciprocalServing'.*ismember(FoodTypeIDs',[14 16]);
%Fats and oils
ConstraintMatrix(7,:)=ReciprocalServing'.*(FoodTypeIDs'==11);
%Sweets
ConstraintMatrix(8,:)=ReciprocalServing'.*(FoodTypeIDs'==25);
%Sodium
ConstraintMatrix(9,:)=FoodMatrix(:,132)'/100;
%Foods not included in DASH diet
ConstraintMatrix(10,:)=(max(ConstraintMatrix(1:8,:))==0);
%-------------------------------------------------------------------------
%Construct the data structure
%-------------------------------------------------------------------------
DASH.ConstraintNames=ConstraintNames;
DASH.ConstraintLBs=ConstraintLBs;
DASH.ConstraintUBs=ConstraintUBs;
DASH.ConstraintMatrix=ConstraintMatrix;
%-------------------------------------------------------------------------
%Finished definition of DASH diet
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
%Diet 9: Vegan (plant-based) diet
%Vegan diet is defined by zero consumption of animal products, eggs and
%dairy products (baked products are also excluded)
%-------------------------------------------------------------------------
ConstraintNames={'Meat, dairy and eggs'};
ConstraintLBs=0;
ConstraintUBs=0;
ConstraintMatrix=double(ismember(FoodTypeIDs',[1 2 3 4 8 9 10 13 15 17 18 19 20 21 22]));
%-------------------------------------------------------------------------
%Construct the data structure
%-------------------------------------------------------------------------
Vegan.ConstraintNames=ConstraintNames;
Vegan.ConstraintLBs=ConstraintLBs;
Vegan.ConstraintUBs=ConstraintUBs;
Vegan.ConstraintMatrix=ConstraintMatrix;
%-------------------------------------------------------------------------
%Finished definition of vegan diet
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
%Diet 10: USDA recommended diet
%-------------------------------------------------------------------------
%USDA recommended diet satisfies all USDA daily nutritional goals
%This diet is age- and gender-dependent. Constraints on this diet include
%lower and upper bounds on nutrient amounts and nutrient ratios
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Define indices of nutrients
%determine the age-gender group based on age_now and gender_now
age_now=29; %My age at the time point of writing this code
gender_now=2; %My gender at the time point of writing this code
ConstraintMatrix=ones(34,8788);

%Create variables
load USDA_Goals.mat
nut_nhanes_pos = [8 9 13 14 21 23 26 31 34 35 37 40 43 44 45 56 65 73 76 77 ...
    82 83 94 102 105 109 111 114 129 132 135 139 140 143 150 157 170 ...
    171 175]; %Indice of NHANES nutrients in USDA SR
%{
Pos_LB = [26 19 9 20 16 17 15 32 2 6 18 37 33 31 24 27 25 14 28 22 38 7 1];
Pos_Sodium=30;
Pos_Fiber=11;
Pos_AMDR=[33 9 13 21 36];
%}
%-------------------------------------------------------------------------
%Initialize part of constraint matrix that is independent on age and gender
%-------------------------------------------------------------------------
%Micronutrients
ConstraintMatrix(1:23,:)=FoodMatrix(:,nut_nhanes_pos(Pos_LB))'/100;
%Calories
ConstraintMatrix(24,:)=FoodMatrix(:,143)'/100;
%Sodium
ConstraintMatrix(25,:)=FoodMatrix(:,nut_nhanes_pos(Pos_Sodium))'/100;
%Fiber/Calories>=0.014
ConstraintMatrix(26,:)=FoodMatrix(:,nut_nhanes_pos(Pos_Fiber))'-0.014*FoodMatrix(:,143)';
%Sugar contribution to calories <=10%
ConstraintMatrix(27,:)=FoodMatrix(:,nut_nhanes_pos(Pos_AMDR(4)))'*4-0.1*FoodMatrix(:,143)';
%Saturated fat contribution to calories <=10%
ConstraintMatrix(28,:)=FoodMatrix(:,nut_nhanes_pos(Pos_AMDR(5)))'*9-0.1*FoodMatrix(:,143)';
%-------------------------------------------------------------------------
%Intialize the AMDR part of constraint matrix that is dependent on age and
%gender
%-------------------------------------------------------------------------
loc_nut_bounds=find(Groups_Gender==gender_now & Groups_AgeUB>=age_now ...
    & Groups_AgeLB<=age_now);
%Lower bound of protein and carbohydrate's contribution to calories
ConstraintMatrix(29:30,:)=FoodMatrix(:,nut_nhanes_pos(Pos_AMDR(1:2)))'*4 ...
    -LBs_AMDR(1:2,loc_nut_bounds)*FoodMatrix(:,143)';
%Lower bound of fat's contribution to calories
ConstraintMatrix(31,:)=FoodMatrix(:,nut_nhanes_pos(Pos_AMDR(3)))'*9 ...
    -LBs_AMDR(3,loc_nut_bounds)*FoodMatrix(:,143)';
%Upper bounds of protein and carbohydrate's contribution to calories
ConstraintMatrix(32:33,:)=FoodMatrix(:,nut_nhanes_pos(Pos_AMDR(1:2)))'*4 ...
    -UBs_AMDR(1:2,loc_nut_bounds)*FoodMatrix(:,143)';
%Upper bound of fat's contribution to calories
ConstraintMatrix(34,:)=FoodMatrix(:,nut_nhanes_pos(Pos_AMDR(3)))'*9 ...
    -UBs_AMDR(3,loc_nut_bounds)*FoodMatrix(:,143)';
%-------------------------------------------------------------------------
%Initialize the upper and lower bounds based on age_now and gender_now
%-------------------------------------------------------------------------
%Micronutrients
ConstraintLBs(1:23)=NutLBs(:,loc_nut_bounds);
ConstraintUBs(1:23)=Inf*ones(23,1);
%Calories
if gender_now == 1
    ConstraintLBs(24)=Cal_Cons_M(Cal_Cons_M(:,1)<=age_now ...
        & Cal_Cons_M(:,2)>=age_now,3);
    ConstraintUBs(24)=Cal_Cons_M(Cal_Cons_M(:,1)<=age_now ...
        & Cal_Cons_M(:,2)>=age_now,4);
else
    ConstraintLBs(24)=Cal_Cons_F(Cal_Cons_F(:,1)<=age_now ...
        & Cal_Cons_F(:,2)>=age_now,3);
    ConstraintUBs(24)=Cal_Cons_F(Cal_Cons_F(:,1)<=age_now ...
        & Cal_Cons_F(:,2)>=age_now,4);
end
%Sodium
ConstraintLBs(25)=0;
ConstraintUBs(25)=Sodium_UB(loc_nut_bounds);
%Bounds for nutrient ratios
ConstraintLBs(26:34)=[0;-Inf;-Inf;0;0;0;-Inf;-Inf;-Inf];
ConstraintUBs(26:34)=[Inf;0;0;Inf;Inf;Inf;0;0;0];
clear age_now gender_now loc_nut_bounds Cal_Cons_F Cal_Cons_M Sodium_UB
%-------------------------------------------------------------------------
%Construct the data structure
%-------------------------------------------------------------------------
USDA.ConstraintNames={'Vitamin E','Calcium','Carbohydrate','Choline','Copper',...
    'Folate','Iron','Magnesium','18:2','18:3','Phosphorus','Potassium',...
    'Protein','Seienium','Vitamin A','Thiamin','Vitamin B-12','Riboflavin',...
    'Vitamin B-6','Vitamin C','Vitamin D','Vitamin K','Zinc',...
    'Calories',...
    'Sodium',...
    'Fiber ratio_lb',...
    'Sugar ratio_ub','Saturated fat ratio_ub',...
    'Protein ratio_lb','Carbohydrate ratio_lb','Fat ratio_lb',...
    'Protein ratio_ub','Carbohydrate ratio_ub','Fat ratio_ub'};
USDA.ConstraintLBs=ConstraintLBs(:);
USDA.ConstraintUBs=ConstraintUBs(:);
USDA.ConstraintMatrix=ConstraintMatrix;
%-------------------------------------------------------------------------
%Finished definition of USDA recommended diet
%-------------------------------------------------------------------------

%Create list of pre-defined diets
DietList={Mediterranean,Japanese,Vegetarian,Vegan,DASH,Paleo,Ketogenic,Atkins,American,USDA};
DietNames={'Mediterranean','Japanese','Vegetarian','Plant-based','DASH','Paleo',...
    'Ketogenic','Atkins','American','USDA'};





