function [expectations theta] = fitGaussMixture(X,k,init,options)
% [expectations theta] = fitGaussMixture(X,k,init,options)
%
%	INPUT:	X: 		n x d data matrix for n data points of dimension d
%			k: 		number of mixture components
%			init: 	method of initialization for the mixture-parameters
%				'k-means': estimate initial params from clustering
%				'Start': use user-defined parameters, located in the struct options.S: .mu,.Sigma,.ComponentProperties

[N dim] = size(X);

% initialization
if strcmp(init,'kmeans')
	Xo = X(sum(isnan(X),2)==0,:);
	[L,C] = kmeans(Xo,k);	
    
	theta.mu = C';
	theta.tau = hist(L,k)./length(L);
	for j = 1:k
		theta.Sigma(:,:,j) = cov(Xo(L==j,:));
	end
	%[theta.mu, theta.Sigma, theta.tau] = kmeansInitMixGauss(X, k);
elseif strcmp(init,'Start')
	theta.mu = options.S.mu;
	theta.tau = options.S.ComponentProportion;
	theta.Sigma = options.S.Sigma;
end

handles.estep = @GMM_EStep;
handles.mstep = @GMM_MStep;

densEst = GMM_densEst(X,theta);
logLike = sum(log(sum(repmat(theta.tau,N,1).*densEst,2)));
fprintf('Log-Likelihood after initialization: %.3f\n',logLike)

[expectations theta] = EMGeneral(X,[],theta,handles,struct());

end
