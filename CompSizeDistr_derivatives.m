function Db_data = CompSizeDistr_derivatives(p,N)
ep = 1e-6; % epsilon for evaluating derivatives using finite differences
plength = length(p);
p = p(1:plength);
p = p/sum(p);
maxdegree = plength-1;
% compute derivative dCC(k)/dp(j) using finite differences  
Mbasis = getONB(plength);
Db = zeros(N,maxdegree);
for j = 1 : maxdegree 
    [~,bplus,~,~,~] = MolSizeDistr_pi(p + ep*Mbasis(:,j),N);
    [~,bminus,~,~,~] = MolSizeDistr_pi(p - ep*Mbasis(:,j),N);
    Db(:,j) = (bplus-bminus)*0.5/ep;
end
Db_data = Db*Mbasis';
end
%%
%%
function Mbasis = getONB(plength)
%% find orthonormal basis in the plain p0 + p1 + ... + p_{maxdegree} = 1
A = [ones(plength,1),[zeros(1,plength-1);eye(plength-1)]];
[Q,~] = qr(A);
Mbasis = Q(:,2:end);
end
