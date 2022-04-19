function [a, b] = posterior_updateBi(y, apre, bpre) 
% update two parameters for beta distribition
    if y>0
        a = apre+1;
        b = bpre;
    else
        a = apre;
        b = bpre+1;
    end