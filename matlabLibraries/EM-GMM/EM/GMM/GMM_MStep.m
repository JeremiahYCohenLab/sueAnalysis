function [theta logLike] = GMM_MStep(X,expectations,thetaOld)
%
% [theta logLike] = GMM_MStep(X,expectations)
%
% Input: 	X - data matrix
%			expectations - posterior of latent indicator variables z
%
% Output: 	theta - struct holding the maximized parameters
%			logLike - value of the log-likelihood function

[N dim] = size(X);
K = size(expectations.z,2);
N_k = sum(expectations.z);

% calculate function value of Q(\theta|\theta^{old})
%Q = GMM_calcQ(X,expectations,thetaOld);
%fprintf('Q before: %.2f',Q);

% update parameters
Xtmp = X;
for j = 1:K
	% handle missing data; replace the missing parts in x by their expectation
	if sum(isnan(X(:))) > 0
		Xtmp(isnan(X)) = expectations.x{j}(isnan(X));
	end
    theta.mu(:,j) = 1/N_k(j)*sum(repmat(expectations.z(:,j),1,dim).*Xtmp);
    theta.Sigma(:,:,j) = 1/N_k(j)*(repmat(expectations.z(:,j),1,dim).*(Xtmp-repmat(theta.mu(:,j)',N,1)))'*(Xtmp-repmat(theta.mu(:,j)',N,1));
    theta.tau(j) = N_k(j)/N;
end

% add A_{nk} = \Lambda_k^{{mm}^{-1}} matrices to Sigma (see script)
if sum(isnan(X(:))) > 0
	pointsMissing = find(sum(isnan(X),2)>0);

	for j = 1:K
		for n = 1:length(pointsMissing)
			componentsMissing = isnan(X(pointsMissing(n),:));
			theta.Sigma(componentsMissing,componentsMissing,j) = theta.Sigma(componentsMissing,componentsMissing,j) + 1/N_k(j)*expectations.z(pointsMissing(n),j)*inv(expectations.xx{j}(componentsMissing,componentsMissing));
		end
	end
end

%Q = GMM_calcQ(X,expectations,theta);
%fprintf(' Q after: %.2f\n',Q);

densEst = GMM_densEst(X,theta);
logLike = sum(log(sum(repmat(theta.tau,N,1).*densEst,2)));

end
