function Q = GMM_calcQ(X,expectations,theta)
%
% calculates the function value for the expectation of the complete-data log-likelihood with respect to the posterior distribution

K = size(expectations.z,2);

Xtmp = X;
Q = 0;
for j = 1:K
    if sum(isnan(X(:))) > 0
        Xtmp(isnan(X)) = expectations.x{j}(isnan(X));
    end
    Q = Q + expectations.z(:,j)'*(log(theta.tau(j))+log(mvnpdf(Xtmp,theta.mu(:,j)',squeeze(theta.Sigma(:,:,j)))));
end

for j = 1:K
	theta.Prec(:,:,j) = inv(squeeze(theta.Sigma(:,:,j)));
end


if sum(isnan(X(:))) > 0
    pointsMissing = find(sum(isnan(X),2)>0);

    for j = 1:K
        for n = 1:length(pointsMissing)
            componentsMissing = isnan(X(pointsMissing(n),:));
			Q = Q + trace(expectations.z(pointsMissing(n),j)*inv(expectations.xx{j}(componentsMissing,componentsMissing))*squeeze(theta.Prec(componentsMissing,componentsMissing,j)));
        end
    end
end


end
