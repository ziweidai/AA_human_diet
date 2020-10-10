%-----------------------------------------------------------------------
% Design MR diets
%-----------------------------------------------------------------------
DefineDiets;
%----------------------------------------------------
%Calculate AA ranges in diets
%----------------------------------------------------
i_met=14; % Index for methionine
i_diet = 10; % Index for USDA recommended diet
n_diets = 100; % Sample 100 MR diets

Diet=DietList{i_diet};
ms=MinServing(CompletePos);
ms(ms<0)=10;

prob_aa_diet.blx=zeros(2335,1);
%No more than 1kg for one food
prob_aa_diet.blc=[Diet.ConstraintLBs;0];
prob_aa_diet.buc=[Diet.ConstraintUBs;3000];
%No more than 3kg total food
prob_aa_diet.a=[Diet.ConstraintMatrix(:,CompletePos);ones(1,length(CompletePos))];
prob_aa_diet.bux=ms*3;

c=FoodMatrix(CompletePos,AAPos(i_met))/100;
d=FoodMatrix(CompletePos,ProteinPos)/100;
prob_aa_diet.c=c;
[~,res_min] = mosekopt('minimize echo(0)',prob_aa_diet);
met_min=res_min.sol.bas.pobjval;
diet_minmet=res_min.sol.bas.xx;
[~,res_max] = mosekopt('maximize echo(0)',prob_aa_diet);
met_max=res_max.sol.bas.pobjval;

diets=zeros(2335,100);
met_levels=zeros(100,1);
prot_levels=zeros(100,1);

prob_mr_diet=prob_aa_diet;
prob_mr_diet.a=[prob_mr_diet.a;c'];
prob_mr_diet.blc=[prob_mr_diet.blc;0];
prob_mr_diet.buc=[prob_mr_diet.buc;0.5];
count=0;
while count<100
    prob_mr_diet.c=c;
    prob_mr_diet.bux=3*ms.*(randi(2,2335,1)-1);
    [~,res] = mosekopt('minimize echo(0)',prob_mr_diet);
    if strcmp(res.sol.bas.solsta,'OPTIMAL')
        count=count+1;
        met_levels(count)=res.sol.bas.pobjval;
        diets(:,count)=res.sol.bas.xx;
        prot_levels(count)=d'*res.sol.bas.xx;
    end
end
[~,idx]=sort(met_levels);
diets=diets(:,idx);
met_levels=met_levels(idx);
prot_levels=prot_levels(idx);

nrows=max(sum(diets~=0));
for i=1:100
    names=repmat({''},nrows,1);
    values=repmat(0,nrows,1);
    n=sum(diets(:,i)~=0);
    names(1:n)=FoodNames(CompletePos(diets(:,i)~=0));
    values(1:n)=nonzeros(diets(:,i));
    if i==1
        t=table(names,values);
        t.Properties.VariableNames=strcat('diet_',num2str(i),{'_foods','_weights'});
    else
        t_add=table(names,values);
        t_add.Properties.VariableNames=strcat('diet_',num2str(i),{'_foods','_weights'});
        t=[t t_add];
    end
end
writetable(t,'MR_Diets.csv');