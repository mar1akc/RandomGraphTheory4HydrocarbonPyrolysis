function [P_distr,pi_distr,uflag,u,S] = MolSizeDistr_pi(p,N)
%% bk = probability for a randomly picked component to be of size k
[P_distr,uflag,u,S] = H0distribution(p,N);
e = (1:N)';
b = (P_distr./e);
Z = sum(b);
pi_distr = b/Z;
end
