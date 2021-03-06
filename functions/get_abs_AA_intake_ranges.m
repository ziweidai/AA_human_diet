%-----------------------------------------------------------------------
% Calculate ranges of AA daily intake in different diets with calories
% level [1800,2200]
%-----------------------------------------------------------------------

Diet_AAmin=zeros(10,18);
Diet_AAmax=zeros(10,18);
label_baby_food = double(contains(FoodNames(CompletePos),'Babyfood'))';
for i=1:10
    Diet=DietList{i};
    DietNames(i)
    prob_aa_diet.blx=zeros(2335,1);
    prob_aa_diet.bux=1000*ones(2335,1); 
    %No more than 1kg for one food
    prob_aa_diet.blc=[Diet.ConstraintLBs;0;-Inf]; 
    prob_aa_diet.buc=[Diet.ConstraintUBs;3000;0];  
    %No more than 3kg total food
    prob_aa_diet.a=[Diet.ConstraintMatrix(:,CompletePos);ones(1,length(CompletePos));label_baby_food];
    if i<9
        prob_aa_diet.blc=[prob_aa_diet.blc;1800];
        prob_aa_diet.buc=[prob_aa_diet.buc;2200];
        prob_aa_diet.a=[prob_aa_diet.a;FoodMatrix(CompletePos,143)'/100];
    end
    for j=1:18
        c=FoodMatrix(CompletePos,AAPos(j))/100;
        prob_aa_diet.c=c;
        [~,res_min] = mosekopt('minimize echo(0)',prob_aa_diet);
        Diet_AAmin(i,j)=res_min.sol.bas.pobjval;
        [~,res_max] = mosekopt('maximize echo(0)',prob_aa_diet);
        Diet_AAmax(i,j)=res_max.sol.bas.pobjval;
    end
end

Diet_AAmin_Norm=Diet_AAmin./max(Diet_AAmin); %Normalize to column maximum
Diet_AAmax_Norm=Diet_AAmax./max(Diet_AAmax);
Diet_AAVar=log10(Diet_AAmax./Diet_AAmin); %Variability of daily AA intake
Diet_AAVar(isinf(Diet_AAVar))=0;
Diet_Var_Pos=find(max(Diet_AAVar')>0);
Diet_AAVar=Diet_AAVar(max(Diet_AAVar')>0,:);

%Plot ranges of absolute AA intake in human dietary patterns
cmap_now = brewermap(100,'RdYlGn');
figure;
colormap(cmap_now(75:-1:26,:));
subplot(1,2,1);
heatmap_cluster(Diet_AAmin_Norm',AANames,DietNames,[0 1]);
title('Minimal daily intake');
subplot(1,2,2);
heatmap_cluster(Diet_AAmax_Norm',AANames,DietNames,[0 1]);
title('Maximal daily intake');
colorbar('Ticks',[0 1],'TickLabels',{'0','Row max'});

figure;
colormap(cmap_now(75:-1:26,:));
heatmap_cluster(Diet_AAVar',AANames,DietNames(Diet_Var_Pos),[0 3]);
title('Variability');
colorbar('Ticks',[0 3]);

clear cmap_now
