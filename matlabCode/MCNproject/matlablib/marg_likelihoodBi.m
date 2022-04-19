function py = marg_likelihoodBi(y, a, b)
% p(y|beta(a,b))
 py = beta(a+y, b+1-y)/beta(a,b);