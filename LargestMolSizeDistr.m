function LargestMolSizeDistr()
Nsamples = 10000; % the number of random graph samples per dataset
fname = 'Data/Degrees_predictions_10reac.csv';
[p_10reac,data_info,NCvector,HCratio,temperature] = read_data(fname);
Ndata = size(p_10reac,1);
largest_mol_hist_data_10RM_RGS = cell(3,Ndata);
for data_index = 1 : Ndata
    p = p_10reac(data_index,:);
    n = NCvector(data_index);
    Mmax = zeros(Nsamples,1);
    for j = 1 : Nsamples
        [G,~] = SampleConfigurationModel(p,n);
        [~,molecules] = DFS(G);
        Mmax(j) = max(molecules);
    end
    Mmax_max = max(Mmax);
    Mmax_min = min(Mmax);
    Mvector = Mmax_min:Mmax_max;
    Mvector_length = length(Mvector);
    Mcounts = zeros(size(Mvector));
    for s = 1 : Mvector_length
        ind = find(Mmax == Mvector(s));
        if ~isempty(ind)
            Mcounts(s) = length(ind);
        end
    end
    largest_mol_hist_data_10RM_RGS{1,data_index} = strcat(data_info(1,data_index),",",num2str(temperature(data_index)),"K,",num2str(n));
    largest_mol_hist_data_10RM_RGS{2,data_index} = Mcounts;
    largest_mol_hist_data_10RM_RGS{3,data_index} = Mvector;
    fprintf('data_index = %d\n',data_index);
end
save('Data/largest_mol_hist_data_10RM_RGS.mat',...
    'largest_mol_hist_data_10RM_RGS','HCratio','temperature');
end