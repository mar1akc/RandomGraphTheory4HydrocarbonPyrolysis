function driver_new_new()
close all
fsz = 20;
% this function reads data from a set of simulations and calls programs
% predicted degree distributions
fname = 'Data_new/Degrees_predictions_10reac.csv';
[p_10reac,data_info] = MakeMatFile4PredDegrees(fname);
NCvector = [256*ones(1,8),250,250,160,288,240,160,320,320,320]';
Ndata = length(NCvector);
p_10reac = [36.58246485	91.39742897	85.19200474	36.23351377	6.594587672
35.67419572	91.44534613	85.50304491	37.01884428	6.358568954
34.81952133	91.05656257	86.41473641	37.59063061	6.118549078
33.86645673	90.96694471	86.89153228	38.33218288	5.942883392
33.00562986	90.65739576	87.41623697	39.13855316	5.78218425
29.52656019	88.84429999	89.87822451	42.27018711	5.480728203
25.59790143	85.46952966	92.63197703	46.7724527	5.528139179
21.86683613	81.90276456	94.64777516	51.86503899	5.717585153
67.66772727	105.330856	59.68225464	15.70194619	1.617215924
54.50852232	104.8284051	68.31514237	20.63590172	1.712028508
88.94639989	56.87191295	12.80572365	1.329856141	0.046107368
27.05658572	87.52357636	104.1068171	56.94862252	12.36439833
64.93346819	102.008626	56.9240093	14.70072109	1.325888621
61.47403042	69.00546623	25.01065778	4.337837059	0.172008517
41.28706794	113.0918192	109.5932238	48.81036936	7.217519669
69.7579815	134.5301278	86.94114018	26.57395208	2.19679845
24.84945644	93.91426531	119.1340334	70.70148947	11.40075534];
p_10reac = p_10reac.*((1./NCvector)*ones(1,5));
% degrees and molecule size distr from MD simulations 
% fname_DegAndMolSize = 'Data/DegreeAndMolSize.mat';
fname_DegAndMolSize = 'Data_new/DegreeAndMolSize.mat';
d = load(fname_DegAndMolSize);
MDorder = [3:10,1:2,12,11,13,17,14,15,16]; % order in which data should be taken
dd = d.DegreeDistribution_cell;
sd = d.mol_size_distr_cell;
% standard deviations
% std_data = load('Data/Degrees_std_equil.mat');
% perr = 2*std_data.ss; % error estimate for p_k is two standard deviations
%% separate data for various kinds of mixtures and various temperatures
% array structure: initial molecules, temperature, experiment #
fid1 = fopen('giant_comp.txt','w');
fid2 = fopen('deg_distr.txt','w');
% gradients of comp size distr w.r.t. degree distr
Dpi_data = CompSizeDistr_derivatives(fname_DegAndMolSize);
%
% data for computing the KL divergence as a function of Temperature and H/C
% ratio
KLdata = zeros(Ndata,5);
for data_index = 1 : Ndata 
    init_mol = cell2mat(data_info(1,data_index));
    if data_index <= 10 || data_index == 12 
    Ncarbon = str2num(init_mol(2));
    Nhydrogen = str2num(init_mol(4:end));
    else
        if data_index == 11 || data_index == 14 % CH4
            Ncarbon = 1;
            Nhydrogen = 4;
        end
        if data_index == 13 % mix
            Ncarbon = 240;
            Nhydrogen = 732;
        end
    end

    temperature = cell2mat(data_info(2,data_index));
    TK = str2num(temperature(1:4));
    p = dd{2,MDorder(data_index)}';
    plength = min(5,find(p>0, 1, 'last' ));
    p = p(1:plength);
    p = p/sum(p);
    maxdegree = plength-1;
    fprintf('Max degree = %d\n',maxdegree);
    % degree distr extracted from MD simulations
    [CC,uflag,u,S] = H0distribution(p(1:plength));
%     fname = strcat('H0distributions/CompSizeDistr',init_mol,temperature,sprintf('_maxdeg%d.mat',maxdegree));
%     data = load(fname);
%     CC = data.CCdistr;
%     uflag = data.uflag;
    if uflag == 1
%         u = data.u;
%         S = data.S;
        fprintf('u = %d, S = %d, sum(CC) = %d\n',u,S,sum(CC));
    else
        u = 1;
        S = 0;
    end
    % degree distr extracted from 10-reaction model
    [CC10,uflag,u10,S10] = H0distribution(p_10reac(data_index,:)');
%     fname = strcat('H0distributions/CompSizeDistr',str1,sprintf('_maxdeg%d.mat',maxdegree));
%     data = load(fname);
%     CC10 = data.CCdistr;
%     uflag = data.uflag;
    if uflag == 1
%         u10 = data.u;
%         S10 = data.S;
        fprintf('u10 = %d, S10 = %d, sum(CC) = %d\n',u10,S10,sum(CC10));
    else
        u10 = 1;
        S10 = 0;
    end
    tname = strcat(init_mol,", ",temperature);
%     l1name = strcat(init_mol,", ",temperature,", RGT");
    l1name = "RGT";
%    l11name = strcat(init_mol,", ",temperature,", MD data");
%    l12name = strcat(init_mol,", ",temperature,", 10RM + RGT");
    l12name = "10RM+RGT";
    % the actual size distribution
    msd = sd{2,MDorder(data_index)};
    msd = msd/sum(msd);
    nC = length(msd);
    msd = msd./(1:nC)';
    msd = msd/sum(msd);
    if uflag == 1
        ngiant = nC*S;
        kmax = min(length(CC),round(nC-ngiant));
        fprintf('<n giant> = %d, kmax = %d\n',round(ngiant),kmax);
    else
        kmax = length(CC);
    end
    klast = find(msd>0, 1, 'last');
    l2name = "MD data";
%     l2name = strcat(init_mol,", ",temperature," MD data");
    pidistr = CC./(1:length(CC))';
    pidistr = pidistr/sum(pidistr);
    pidistr10 = CC10./(1:length(CC10))';
    pidistr10 = pidistr10/sum(pidistr10);
    aux = nC*cumsum(CC);
    if uflag == 1
        pilast = find(aux<kmax,1,'last');
    else
        inonzero = find(msd>0,1,'last'); 
        pilast = min(kmax,inonzero);
    end
    if data_index == 15
        pilast = 20;
    end
    fprintf('pilast = %d\n',pilast);
    % error bars for predicted distribution
    Dpi = squeeze(Dpi_data(MDorder(data_index),:,:)); 
    std4deg = dd{3,MDorder(data_index)}';
    perr = 2*std4deg(1:5);
     ErrBarPos = abs(Dpi)*perr;
    ErrBarNeg = min(ErrBarPos,pidistr*0.9999);
    % error bars for distribution extracted from md data
    pilast = min(20,pilast); %length(CC);
    msd_std = sd{3,MDorder(data_index)};
    msd_err_pos = 2*msd_std(1:pilast)/(NCvector(data_index)*sum(pidistr));
    msd_err_neg = min(msd_err_pos,pidistr(1:pilast)*0.9999);
%     msd(1:pilast);
%     size(msd(1:pilast))
%     size(msd_err_pos)
%     size(msd_err_neg)
    % plotting
    fig = ceil(data_index+1);
    figure(fig); 
    hold on
%     for k = 1 : pilast
%         fprintf('%d\t%d\t%d\t%d\n',k,msd(k),msd_err_neg(k),msd_err_pos(k));
%     end
    errorbar(1:pilast,msd(1:pilast),msd_err_neg,msd_err_pos,...
        '.','Markersize',30,'color','k','Displayname',l2name);
%    plot(msd,...
%        '.','Markersize',30,'color','k','Displayname',l2name);
    errorbar(1:pilast,pidistr(1:pilast),...
        ErrBarNeg(1:pilast),ErrBarPos(1:pilast),...
        '--','Linewidth',1,'color','b','Displayname',l1name);
    plot(1:pilast,pidistr10(1:pilast),...
        'Linewidth',1,'color','r','Displayname',l12name);
%    if mod(data_index,2) == 0 || data_index == 23
        legend
        xlabel('k','fontsize',fsz);
        ylabel('\pi_k','fontsize',fsz);
        title(tname);
        axis([1,pilast,msd(pilast)*0.9,1])
%         if data_index == 23
%              axis([1,10,1e-7,1]);
%         end
        set(gca,'YScale','log','Xscale','log','fontsize',fsz);
        grid
        figname = strcat('Figures/',init_mol,'_',temperature,'_new.eps');
        saveas(gcf,figname,'epsc');
%    end 
    % latex table
    fprintf(fid1,'%d & %s & %s & %d & %.4f & %.4f & %d\n',data_index,init_mol,temperature,nC,u,S,round(nC*S));
    fprintf(fid2,'%d & %s & %s & %.4f & %.4f & %.4f & %.4f & %.4f\n',data_index,init_mol,temperature,p(1),p(2),p(3),p(4),p(5));
    fprintf('Data index = %d\n',data_index)
    Ncarbon
    Nhydrogen
    TK
    init_mol
%     pi_data = msd(1:pilast);
%     pi_rgt = pidistr(1:pilast);
%     pi_10rmgrt = pidistr10(1:pilast);
%         KL(msd(1:pilast),pidistr(1:pilast))
%         KL(msd(1:pilast),pidistr10(1:pilast))
    KLdata(data_index,:) = [Ncarbon,Nhydrogen,TK,...
        W1(msd(1:pilast),pidistr(1:pilast)),...
        W1(msd(1:pilast),pidistr10(1:pilast))];
end
fclose(fid1);
fclose(fid2);
save('KLdata_new.mat','KLdata');
end
%%
function [CCdistr,uflag,u,S] = H0distribution(p)
tol = 1e-14;
plength = length(p);
maxdegree = plength-1;
kvector = (0:maxdegree)';
q = zeros(maxdegree,1); % excess degree distribution generated by G1
meandegree = kvector'*p; % mean degree
q = kvector(2:end).*p(2:end)/meandegree; % qk = (k+1)p_{k+1}/meandegree
%
N = 128; % the max molecule size
%% compute H0 depending on max degree
uflag = 0;
u = 1;
S = 0;
if maxdegree == 2 % chain model, completely analytic   
    % compute derivatives analytically
    C = zeros(N,1); % component size distribution
    C(1) = p(1); %p0;
    C(2) = p(2)*q(1);%p1*q0;
    C(3) = p(2)*q(1)*q(2) + p(3)*q(1)^2; %p1*q0*q1+p2*q0^2;
    C(4:N) = q(1)*q(2).^(1:N-3).*(p(2)*q(2)+(2:N-2)*p(3)*q(1));   %q0*q1.^(1:N-3).*(p1*q1+(2:N-2)*p2*q0);
else
    if maxdegree == 3
        G0 = @(x)p(1)+p(2)*x+p(3)*x.^2+p(4).*x.^3;
        G1 = @(x) q(1) + x*q(2) + x.^2*q(3);
        H1 = @(x)2*x*q(1)./((1-x*q(2))+sqrt((1-x*q(2)).^2-4*x.^2*q(1)*q(3)));
        H0 = @(x)x.*G0(H1(x));
        % check if there is a giant component
        giantfun = -p(2) + 3*p(4);
        if giantfun < 0
            fprintf('There is no giant component\n');
        else
            root = sort(roots([q(3),q(2)-1,q(1)]),'ascend');
            u = root(1);
            fprintf('There is a giant component: u = %d\n',u);
            S = 1-G0(u);
            fprintf('Fraction in the giant component: %d\n',S);
            uflag = 1;
        end
    else % maxdegree = 4
        G0 = @(x)p(1) + p(2)*x + p(3)*x.^2 + p(4)*x.^3 + p(5)*x.^4;
        G1 = @(x)q(1) + q(2)*x + q(3)*x.^2 + q(4)*x.^3;
%        xG0 = @(x)x.*p(1) + p(2)*x.^2 + p(3)*x.^3 + p(4)*x.^4 + p(5)*x.^5;
        H1 = @(x)cubicroot(x,q);
        H0 = @(x)x.*G0(H1(x));
        % check if there is a giant component
        giantfun = -p(2) + 3*p(4) + 8*p(5);
        if giantfun < 0
            fprintf('There is no giant component\n');
        else
            root = sort(roots([q(4),q(3),q(2)-1,q(1)]),'ascend');
            ind = find(root>0 & root<1+tol);
            u = root(ind(1));
            fprintf('There is a giant component: u = %d\n',u);
            S = 1-G0(u);
            fprintf('Fraction in the giant component: %d\n',S);
            fprintf('u = H1(1) = %d, G0(u) = %d, 1-G0(u) = %d,  H0(1) = %d, 1-H0(1) = %d\n',H1(1),G0(u),1-G0(u), H0(1),1-H0(1));
            fprintf('H0(1) = x*G0(H1(x) = %d for x = 1\n',G0(H1(1)));
            uflag = 1;
        end
    end
     % compute derivatives of H0 using Cauchy formula
    npt = 1000;
    t = linspace(0,2*pi,npt+1)';
    tm = linspace(0,2*pi,2*npt+1)';
    tm(1:2:end) = [];
    r = 1;
    z = r*exp(1i*t);
    dz = circshift(z,[-1,0]) - z;
    dz(end) = [];
    zm = r*exp(1i*tm);
    C = zeros(N,1);
    for n = 1 : N
        C(n) = sum(dz.*H0(zm)./zm.^(n+1))*0.5/(pi*1i);
    end
end
CCdistr = real(C);
% fname = strcat('H0distributions/CompSizeDistr',str,sprintf('_maxdeg%d.mat',maxdegree));
% save(fname,'CCdistr','uflag','u','S');
end

%%
function root = cubicroot(x,q)
a = q(4)*x;
b = q(3)*x;
c = q(2)*x-1;
d = q(1)*x;
D0 = b.^2-3*a.*c;
D1 = 2*b.^3-9*a.*b.*c + 27*a.^2.*d;
C = (0.5*(D1+sqrt(D1.^2-4*D0.^3))).^(1/3);
xi = -0.5 + 1i*0.5*sqrt(3);
xi2 = xi^2;
r0 = -(b+C+D0./C)./(3*a);
r1 = -(b+xi*C+D0./(C*xi))./(3*a);
r2 = -(b+C*xi2+D0./(C*xi2))./(3*a);
r012 = [r0,r1,r2];
% size(r012)
[~,jmin] = min(abs(r012),[],2);
imin = (1:length(x))';
ind = sub2ind([length(x),3],imin,jmin);
root = r012(ind);
end
%%
%%
function d = KL(p,q)
ind_p = find(p>0);
ind_q = find(q>0);
ind = intersect(ind_p,ind_q);
p(ind) = p(ind)/sum(p(ind));
q(ind) = q(ind)/sum(q(ind));
d = sum(p(ind).*log(p(ind)./q(ind)));
end
%%
function d = W1(p,q)
p = p/sum(p);
q = q/sum(q);
pcum = cumsum(p);
qcum = cumsum(q);
d = sum(abs(pcum-qcum));
end






