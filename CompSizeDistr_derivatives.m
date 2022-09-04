function Db_data = CompSizeDistr_derivatives(fname)
close all
fsz = 20;
% this function reads data from a set of simulations and calls programs
% computing degree distributions
d = load(fname);
dd = d.DegreeDistribution_cell;
sd = d.mol_size_distr_cell;
sz = size(dd);
Ndata = sz(2);
ep = 1e-6; % epsilon for evaluating derivatives using finite differences

%% separate data for various kinds of mixtures and various temperatures
% array structure: initial molecules, temperature, experiment #
data_info = cell(3,Ndata);
N = 128; % the largest molecule size that is not in the giant component
for data_index = 1 : Ndata
    ss = strsplit(dd{1,data_index},'_');
    for j = 1 : length(ss)
        data_info{j,data_index} = ss{j};
    end    
end
col = lines(Ndata);
% fid1 = fopen('giant_comp.txt','w');
% fid2 = fopen('deg_distr.txt','w');
Db_data = zeros(Ndata,N,5);
figure;
for data_index = 1 : Ndata 
    init_mol = cell2mat(data_info(1,data_index));
    temperature = cell2mat(data_info(2,data_index));
    run = cell2mat(data_info(3,data_index));
    % string indicating chemical composition and temperature
    str = dd{1,data_index};
    % degree distribution of length 6
    p = dd{2,data_index};
    p = p';
    plength = min(5,find(p>0, 1, 'last' ));
    p = p(1:plength);
    p = p/sum(p);
    maxdegree = plength-1;
    fprintf('%s, T = %s, Max degree = %d\n',init_mol,temperature,maxdegree);
    % compute derivative dCC(k)/dp(j) using finite differences  
    Mbasis = getONB(plength);
    [CC,b,u,S] = Distr_pi(p,N);
    fprintf('u = %d, S = %d, sum(CC) = %d\n',u,S,sum(CC));
    Db = zeros(N,maxdegree);
    DS = zeros(1,maxdegree);
    for j = 1 : maxdegree 
         [~,bplus,~,Splus] = Distr_pi(p + ep*Mbasis(:,j),N);
        [~,bminus,~,Sminus] = Distr_pi(p - ep*Mbasis(:,j),N);
        Db(:,j) = (bplus-bminus)*0.5/ep;
        DS(j) = (Splus-Sminus)*0.5/ep;
    end
    Dbnorm = sqrt(sum(Db.^2,2));
    DSnorm = norm(DS);
    if S > 0 
        fprintf('Norm(DS) = %d\n',DSnorm);
    end
    Dbdir = Db.*((1./Dbnorm)*ones(1,maxdegree))*Mbasis';
    Db_data(data_index,:,:) = Db*Mbasis';
%     for j = 1 : 10
%         for k = 1 : plength
%             fprintf('%d\t',Dbdir(j,k));
%         end
%         fprintf('; Dbnorm = %d\n',Dbnorm(j));
%     end
    fig = ceil(data_index);
    subplot(4,5,data_index); 
    hold on
    l1name = strcat(init_mol,", ",temperature);
    % plotting
    plot(Dbnorm./b,'.','Markersize',20,'Displayname',l1name);
        legend('Location','southeast');
        xlabel('k','fontsize',fsz);
        ylabel('||\nabla \pi(k)||/\pi_k','fontsize',fsz);
        set(gca,'fontsize',fsz);
%        set(gca,'YScale','log','Xscale',xscale1,'fontsize',fsz);
        grid
%         figname = strcat('Figures/',init_mol,'_',temperature,'.eps');
%         saveas(gcf,figname,'epsc');
%     end 
    % latex table
%     fprintf(fid1,'%d & %s & %s & %d & %.4f & %.4f & %d\n',data_index,init_mol,temperature,nC,u,S,round(nC*S));
%     fprintf(fid2,'%d & %s & %s & %.4f & %.4f & %.4f & %.4f & %.4f\n',data_index,init_mol,temperature,p(1),p(2),p(3),p(4),p(5));
end
% fclose(fid1);
% fclose(fid2);
%save('pi_grad.mat','Db_data');

end

%%
function Mbasis = getONB(plength)
%% find orthonormal basis in the plain p0 + p1 + ... + p_{maxdegree} = 1
A = [ones(plength,1),[zeros(1,plength-1);eye(plength-1)]];
[Q,~] = qr(A);
Mbasis = Q(:,2:end);
end
