function W1dist4SmallMolSizeDistr()
close all
fsz = 20;
data = load('W1data.mat');
W1 = data.W1data;
Ndata = 14;
HCratio = W1(:,1);
temperature = W1(:,2);
W1_data_rgt = W1(:,3);
W1_data_10rmrgt = W1(:,4);
%%
T = 3200:100:5000;
nT = length(T);
col = jet(nT);
figure(1); hold on; grid;
for j = 1 : 8
    jcol = temperature(j)/100 - 31;
    plot(HCratio(j),W1_data_rgt(j),'.','Markersize',60,'color',col(jcol,:));
end
for j = 9 : 14
    jcol = temperature(j)/100 - 31;
    plot(HCratio(j),W1_data_rgt(j),'p','Markersize',20,...
        'MarkerEdgeColor',col(jcol,:),'MarkerFaceColor',col(jcol,:));
end
axis([2,4,0,1])
colorbar
colormap jet
caxis([3200,5000]);
xlabel('H/C ratio','FontSize',fsz);
ylabel('Wasserstein W_1 distance','FontSize',fsz);
set(gca,'FontSize',fsz);
title('MD data, RGT','FontSize',fsz)
%
figure(2); hold on; grid;
for j = 1 : 8
    jcol = temperature(j)/100 - 31;
    plot(HCratio(j),W1_data_10rmrgt(j),'.','Markersize',60,'color',col(jcol,:));
end
for j = 9 : 14
    jcol = temperature(j)/100 - 31;
    plot(HCratio(j),W1_data_10rmrgt(j),'p','Markersize',20,...
        'MarkerEdgeColor',col(jcol,:),'MarkerFaceColor',col(jcol,:));
end
axis([2,4,0,1])
colorbar
colormap jet
caxis([3200,5000]);
xlabel('H/C ratio','FontSize',fsz);
ylabel('Wasserstein W_1 distance','FontSize',fsz);
set(gca,'Fontsize',fsz)
title('MD data, 10RM + RGT','FontSize',fsz)
end
