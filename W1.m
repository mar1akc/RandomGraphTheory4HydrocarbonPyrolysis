function d = W1(p,q)
p = p/sum(p);
q = q/sum(q);
pcum = cumsum(p);
qcum = cumsum(q);
d = sum(abs(pcum-qcum));
end
