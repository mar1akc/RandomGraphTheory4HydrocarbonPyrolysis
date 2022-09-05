function [G,deg_distr] = SampleConfigurationModel(p,n)
% Input:
% p = the desired degree distribution
% n = the number of vertices
% Output:
% G = adjacency matrix of the random graph
%% Assign degrees to vertices according to p
K = length(p);
pcounts = n*p;
pcount = zeros(size(pcounts));
for j = 2 : K
    pfloor = floor(pcounts(j));
    eta = rand;
    if pfloor+eta < pcounts(j)
        pcount(j) = pfloor;
    else
        pcount(j) = pfloor + 1;
    end
end
pcount(1) = n - sum(pcount(2:K));
aux = cumsum(pcount);
% p(1) = probability of degree 0
w = zeros(n,1);
prev = 1;
for j = 1 : K
    if aux(j) >= prev
        w(prev:aux(j)) = j-1;
        prev = aux(j)+1;
    end
end
% w = vector of degrees
wcumsum = cumsum(w);
m = wcumsum(end); % (the total number of edges)*2
% m must be even. Make it even if it is odd
if mod(m,2) == 1
    w(aux(1)) = w(aux(1)) + 1;
    wcumsum = cumsum(w);
end
m = wcumsum(end);
%% Enumerate stubs and assign the to vertices according to their degrees
stubs = zeros(m,1);
prev = 1;
for j = 1 : n
    stubs(prev:wcumsum(j)) = j-1;
    prev = wcumsum(j) + 1;
end
%% Randomly match stubs
s0 = randperm(m);
s1 = s0(1:m/2);
s2 = s0(m/2+1:end);
edges = zeros(m/2,2);
edges(:,1) = stubs(s1);
edges(:,2) = stubs(s2);
% remove selfloops
ind = edges(:,1)==edges(:,2);
edges(ind,:) = [];
% remove repeated edges
edges = sort(edges,2);
edges = unique(edges,'rows');
% form the adjacency matrix
G = sparse(edges(:,1),edges(:,2),ones(size(edges,1),1),n,n);
G = spones(G + G');
% test degree distribution
deg = sum(G,2);
ndeg = zeros(K,1);
for j = 1 : K
    ndeg(j) = length(find(deg == j-1));
end
sdeg = sum(ndeg);
deg_distr = ndeg/sdeg;
end