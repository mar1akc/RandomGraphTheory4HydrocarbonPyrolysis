function LargestMolHistograms()
close all
data_MD = load("Data/largest_mol_hist_data_MD.mat");
data_10RM_RGT = load("Data/largest_mol_hist_data_10RM_RGT.mat");

d1 = data_MD.largest_mol_hist_data_MD;
d2 = data_10RM_RGT.largest_mol_hist_data_10RM_RGT;
HCratio = data_10RM_RGT.HCratio;
temperature = data_10RM_RGT.temperature;

Ndata = size(d1,2);

%% plot histograms
for i = 1 : Ndata

    figure()
    bar([d1{3,i}]-0.5,[d1{2,i}]/sum([d1{2,i}]),'BarWidth',1,'FaceAlpha',1,'EdgeColor','none')
    hold on
    bar([d2{3,i}]-0.5,[d2{2,i}]/sum([d2{2,i}]),'BarWidth',1,'FaceAlpha',0.5,'EdgeColor','none')

    title(d2{1,i})

    xlabel('Molecule Size: # C-atoms, s')
    ylabel('\zeta_{s}')
    legend('MD Data','10RM + RGS')
    set(gca,'FontSize',20)

    print(strcat('LargestMolHist',d2{1,i}),'-dpdf')

end

%% compute means and standard deviations
for i = 1 : Ndata
    p1 = [d1{2,i}]/sum([d1{2,i}]);
    LMmean = sum(p1.*[d1{3,i}]);
    LMstd = sqrt(sum(p1.*([d1{3,i}]-LMmean).^2));
    fprintf(strcat('MD:       ',d2{1,i},': '));
    fprintf('mean = %.2f, std = %.2f\n',LMmean,LMstd);
    p2 = [d2{2,i}]/sum([d2{2,i}]);
    LMmean = sum(p2.*[d2{3,i}]);
    LMstd = sqrt(sum(p2.*([d2{3,i}]-LMmean).^2));
    fprintf(strcat('10RM+RGT: ',d2{1,i},': '));
    fprintf('mean = %.2f, std = %.2f\n\n',LMmean,LMstd);
end
end