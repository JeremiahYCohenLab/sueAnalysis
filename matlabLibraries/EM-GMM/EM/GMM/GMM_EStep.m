function expectations = GMM_EStep(X,theta)
%
% expectations = GMM_EStep(X,theta)
%
% Input:    X - data matrix
%           theta - struct holding the current parameters
%
% Output:   expectations - expectations for the latent variables

[N dim] = size(X);
K = length(theta.tau);
   
densEst = GMM_densEst(X,theta);
pointsMissing = find(sum(isnan(X),2)>0);
numElemsMissing = sum(sum(isnan(X)));

% update posterior probabilities p(z|x,\beta)
for j = 1:K
	expectations.z(:,j) = theta.tau(j)*densEst(:,j)./sum(repmat(theta.tau,N,1).*densEst,2);
end

% for the subsequent calculations, we need precision matrices
if length(pointsMissing) > 0
	for j = 1:K
		theta.Prec(:,:,j) = inv(squeeze(theta.Sigma(:,:,j)));
	end
end

% calculate the second factor of E[z_{nk}x_n^m] = E[z_{nk}]E[x_n^m]
for j = 1:K
	I = zeros(numElemsMissing,1);
	J = zeros(numElemsMissing,1);
	S = zeros(numElemsMissing,1);

	counter = 1;
	for n = 1:length(pointsMissing)
		missingComponents = isnan(X(pointsMissing(n),:));
		mat = inv(squeeze(theta.Prec(missingComponents,missingComponents,j)))*squeeze(theta.Prec(missingComponents,~missingComponents,j));
		S(counter:counter+sum(missingComponents)-1) = theta.mu(missingComponents,j) - mat*(X(pointsMissing(n),~missingComponents)-theta.mu(~missingComponents,j));
		I(counter:counter+sum(missingComponents)-1) = pointsMissing(n);
		J(counter:counter+sum(missingComponents)-1) = find(missingComponents);
		counter = counter+sum(missingComponents);
	end
	expectations.x{j} = sparse(I,J,S,N,dim);
end

% do not calculate the second factor of E[z_{nk}x_n^m x_n^{mT}] = E[z_{nk}]E[x_n^m x_n^{mT}]
% instead calculate it as required in the M-Step; but precalculate the precision matrices based on the current parameters
for j = 1:K
	expectations.xx{j} = inv(squeeze(theta.Sigma(:,:,j)));
end

end % end function
