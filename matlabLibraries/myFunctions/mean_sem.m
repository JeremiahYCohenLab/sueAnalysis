function [MEAN SEM]=mean_sem(r,dim);
% r: the data mat
% caculate SEM MEAN
% SEM=std/sqrt(n)
% n: data number
% delete NaN data to caculate
% dim=1: work as column
% dim=2: work as row
% by Darry 2019-4-23
[m n]=size(r);
if dim==1% column
    for i=1:n
        x=r(:,i);
        TF = ~isnan(x);
        x=x(TF);
        L=length(x);
        SEM(i)=std(x)/sqrt(L);
        MEAN(i)=mean(x);
    end
elseif dim==2% row
    for i=1:m
        x=r(i,:);
        TF = ~isnan(x);
        x=x(TF);
        L=length(x);
        SEM(i)=std(x)/sqrt(L);
        MEAN(i)=mean(x);
    end
end
clear m n i x L TF 