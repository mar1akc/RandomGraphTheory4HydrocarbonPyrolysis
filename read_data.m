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
