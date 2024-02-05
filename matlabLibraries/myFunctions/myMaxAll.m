function [m, ind] = myMaxAll(A)
% returns max value position as subscription
% A is a D dimentional matrix
% m is maximum value
% ind an 1XD array, indicating position of maximum value in A

[m, I] = max(A,[],'all', 'linear');
sizeA = size(A);
ind = zeros(size(sizeA));
numA = prod(sizeA);

i = length(sizeA) - 1;
remain = I;
while i > 0
    numA = numA/sizeA(i+1);
    quot = ceil(remain/numA);
    remain = remain - (quot-1)*numA;
    ind(i+1) = quot;
    if i == 1
        ind(i) = remain;
    end
    i = i-1;   
end

