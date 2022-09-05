function [N_con_comp,components] = DFS(G)
n = size(G,1);
% white = 0, gray = 1, black = 2
colors = zeros(n,1);
parent = -1*ones(n,1);
d = zeros(n,1);
f = zeros(n,1);
t = 0;
N_con_comp = 0;
for j = 1 : n
    if colors(j) == 0
        N_con_comp = N_con_comp + 1;
        black = length(find(colors == 2));
        [colors,parent,d,f,t] = DFSvisit(G,colors,parent,d,f,j,t);
        newblack = length(find(colors == 2));
        components(N_con_comp) = newblack - black;
    end
end
end
%%
function [colors,parent,d,f,t] = DFSvisit(G,colors,parent,d,f,j,t)
t = t + 1;
d(j) = t;
colors(j) = 1;
ind = find(G(j,:) == 1);
if ~isempty(ind)
    nj = length(ind);
    for k = 1 : nj
        v = ind(k);
        if colors(v) == 0
            parent(v) = j;
            [colors,parent,d,f,t] = DFSvisit(G,colors,parent,d,f,v,t);
        end
    end
end
colors(j) = 2;
t = t + 1;
f(j) = t;
end
     