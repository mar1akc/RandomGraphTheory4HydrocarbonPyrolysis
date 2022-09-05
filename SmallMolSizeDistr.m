function SmallMolSizeDistr()
close all
fsz = 20; % font size for graphics
%% Computes small molecule size distribution using degree distribution
% extracted from MD data and predicted by ten-reaction model and compares
% it with the molecule size distribution extracted from MD data
NCmax = 128; % maximal molecule size 
kmax = 20; % truncation for displaying small molecule size distribution
%% Read data from files
% degree distribution from ten-reaction model
fname = 'Data/Degrees_predictions_10reac.csv';
[p_10reac,data_info,NCvector,HCratio,temperature] = read_data(fname);
Ndata = size(p_10reac,1); % the number of datasets
% degrees and molecule size distr from MD simulations 
fname_DegAndMolSize = 'Data/DegreeAndMolSizeMDdata.mat';
d = load(fname_DegAndMolSize);
dd = d.DegreeDistribution_cell;
% dd = 3-by-17 cell array
% dd{1,j} = string with the name of the dataset j
% dd{2,j} = degree distribution for dataset j extracted from MD simulations
% dd{3,j} = standard deviations for degree distribution j extracted from MD
% simulations
sd = d.MolSizeDistribution_cell;
% sd = 3-by-17 cell array
% sd{1,j} = string with the name of the dataset j
% sd{2,j} = molecule size distribution for dataset j extracted from MD simulations
% sd{3,j} = standard deviations for molecule counts for dataset j 
% set up data array for computing W1 distances 
W1data = zeros(Ndata,4);
% loop over available datasets
for data_index = 1 : Ndata 
    fprintf('Data index = %d\n',data_index)
    init_mol = cell2mat(data_info(1,data_index));
    % p = degrees extracted from MD simulations
    p = dd{2,data_index}';
    plength = min(5,find(p>0, 1, 'last' ));
    p = p(1:plength);
    p = p/sum(p);
    maxdegree = plength-1;
    nC = NCvector(data_index);
    fprintf('Max degree = %d\n',maxdegree);
    % compute molecule size distribution for degree distribution extracted
    % from MD simulations
    [P_distr_RGT,pi_distr_RGT,uflag,u,S] = MolSizeDistr_pi(p(1:plength),NCmax);
    if uflag == 1 % there is a giant component
        % S = fraction in giant component
        % u = the minimal positive solution to u = G1(u)
        fprintf('Deg. Distr. from MD: u = %d, S = %d, sum(P_distr_RGT) = %d\n',u,S,sum(P_distr_RGT));
    else % no giant component
        u = 1;
        S = 0;
    end
    % compute molecule size distribution for degree distribution predicted
    % by ten-reaction model
    [P_distr_10RM_RGT,pi_distr_10RM_RGT,uflag,u10,S10] = MolSizeDistr_pi(p_10reac(data_index,:)',NCmax);
    if uflag == 1
        fprintf('DEg. distr. by 10RM: u = %d, S = %d, sum(P_distr_10RM_RGT) = %d\n',u10,S10,sum(P_distr_10RM_RGT));
    else
        u10 = 1;
        S10 = 0;
    end
    tname = strcat(init_mol,", ",num2str(temperature(data_index)),"K"); % figure title
    l1name = "RGT"; % legend name for RGT
    l2name = "10RM+RGT"; % legend name for 10RM+RGT
    l3name = "MD data"; % legend name for MD data
    % the molecule size distribution extracted from MD simulations
    pi_distr_MD = sd{2,data_index};
    pi_distr_MD = pi_distr_MD/sum(pi_distr_MD);
    if uflag == 1
        ngiant = NCvector(data_index)*S;
        fprintf('<n giant> = %d\n',round(ngiant));
    end
    % truncate distribution at k = kmax
    inonzero = find(pi_distr_MD>0,1,'last'); 
    pi_last = min(kmax,inonzero);
    fprintf('pi_last = %d\n',pi_last);
    % error bars for predicted molecule size distribution for degree 
    % distribution extracted from MD data
    Dpi = CompSizeDistr_derivatives(p,NCmax);
    std4deg = dd{3,data_index}'; % standard deviations for degrees
    perr = 2*std4deg(1:5); % error = 2*std
    ErrBarPos = abs(Dpi)*perr;
    ErrBarNeg = min(ErrBarPos,pi_distr_RGT*0.9999);
    % error bars for distribution extracted from md data
    pi_distr_MD_std = sd{3,data_index};
    pi_distr_MD_err_pos = 2*pi_distr_MD_std(1:pi_last)/(NCvector(data_index)*sum(pi_distr_RGT));
    pi_distr_MD_err_neg = min(pi_distr_MD_err_pos,pi_distr_MD(1:pi_last)*0.9999);
    % plotting
    figure(data_index); 
    hold on
    errorbar(1:pi_last,pi_distr_MD(1:pi_last),...
        pi_distr_MD_err_neg,pi_distr_MD_err_pos,...
        '.','Markersize',30,'color','k','Displayname',l3name);
    errorbar(1:pi_last,pi_distr_RGT(1:pi_last),...
        ErrBarNeg(1:pi_last),ErrBarPos(1:pi_last),...
        '--','Linewidth',1,'color','b','Displayname',l1name);
    plot(1:pi_last,pi_distr_10RM_RGT(1:pi_last),...
        'Linewidth',1,'color','r','Displayname',l2name);
    legend
    xlabel('s','fontsize',fsz);
    ylabel('\pi_s','fontsize',fsz);
    title(tname);
    axis([1,pi_last,pi_distr_MD(pi_last)*0.9,1])
    set(gca,'YScale','log','Xscale','log','fontsize',fsz);
    grid
    figname = strcat(tname,'.eps');
    saveas(gcf,figname,'epsc');
    fprintf('Data index = %d\n',data_index)
    W1data(data_index,:) = [HCratio(data_index),temperature(data_index),...
        W1(pi_distr_MD(1:pi_last),pi_distr_RGT(1:pi_last)),...
        W1(pi_distr_MD(1:pi_last),pi_distr_10RM_RGT(1:pi_last))];
end
save('W1data.mat','W1data');
end
%%

%%
function d = W1(p,q)
p = p/sum(p);
q = q/sum(q);
pcum = cumsum(p);
qcum = cumsum(q);
d = sum(abs(pcum-qcum));
end

%%
function [deg_distr,data_info,NCvector,HCratio,temperature] = read_data(fname)

data = readtable(fname);
txt = data.Var1;

d0 = data.Degree0;
d1 = data.Degree1;
d2 = data.Degree2;
d3 = data.Degree3;
d4 = data.Degree4;

NCvector = data.NC;
HCratio = data.HCratio;
temperature = data.Temperature;

Ndata = size(d0,1);
data_info = cell(2,Ndata);
deg_distr = zeros(Ndata,5);
j = 0;
count = 0;
while j < Ndata
    j = j + 1;
    count = count + 1;
    t = txt{j};
    tsplit = strsplit(t,{'_'});
    t1 = tsplit{1};
    t2 = tsplit{2};
    
    p = [d0(j),d1(j),d2(j),d3(j),d4(j)];
    p = p/sum(p);
    deg_distr(count,:) = p;
    data_info{1,count} = t1;
    data_info{2,count} = t2;
end
deg_distr(count+1:end,:) = [];
end



