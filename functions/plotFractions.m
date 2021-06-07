function plotFractions(data,samp_names,component_names)
%Draw rectangles to visualize fractions of different components in a whole
%Input arguments:
% data - n_s*n_c matrix, each row is a sample, each column a component,
% sum(row) = 1
% samp_names - cell array storing names of samples
% component_names - cell array storing names of components
[n_s,n_c]=size(data);
%{
colors=[brewermap(8,'Accent');brewermap(8,'Dark2');brewermap(12,'Paired');...
    brewermap(9,'Pastel1');brewermap(8,'Pastel2');brewermap(9,'Set1');...
    brewermap(8,'Set2');brewermap(12,'Set3')];
%}
colors=[brewermap(9,'Pastel1');brewermap(8,'Pastel2');...
    brewermap(8,'Set2');brewermap(12,'Set3')];
rng(2);
rp=randperm(size(colors,1));
colors=colors(rp,:);

dy=0.1;ygap=0.05;x0=0;y0=0;y1=y0+dy;
%figure;
hold on;
ytickpos=zeros(1,n_s);
D=pdist(data);
tree=linkage(D,'average');
leafOrder_row=optimalleaforder(tree,D);
%leafOrder_row=1:n_s;
for i=1:n_s
    ytickpos(i)=(y0+y1)/2;
    for j=1:n_c
        x1=x0+data(leafOrder_row(i),j);
        h=fill([x0 x1 x1 x0 x0],[y0 y0 y1 y1 y0],colors(j,:));
        %set(h,'EdgeColor','none');
        x0=x1;
    end
    x0=0;
    y0=y1+ygap;
    y1=y0+dy;
end
xtickpos=zeros(1,n_c);
if n_s>100
    for i=1:n_c
        xtickpos(i)=sum(mean(data(:,1:i)))-0.5*mean(data(:,i));
    end
else
    for i=1:n_c
        xtickpos(i)=sum((data(leafOrder_row(1),1:i)))...
            -0.5*(data(leafOrder_row(1),i));
    end
end
xticks(xtickpos);
xticklabels(component_names);
if length(samp_names)>0
    yticks(ytickpos);
    yticklabels(samp_names(leafOrder_row));
else
    yticks([]);
    yticklabels([]);
end
xtickangle(45);
xlim([0 1]);
ylim([0 y0-ygap]);
drawnow;
hAxes=gca;
hAxes.XRuler.Axle.LineStyle = 'none';
hAxes.YRuler.Axle.LineStyle = 'none';
end