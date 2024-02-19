function densEst = GMM_densEst(X,theta)

[N dim] = size(X);
K = length(theta.tau);

% data points that are fully observed
pointsObserved = sum(isnan(X),2)==0;
densEst = zeros(N,K);

% estimate density for observed data points
for j = 1:K
    densEst(pointsObserved,j) = mvnpdf(X(pointsObserved,:),theta.mu(:,j)',squeeze(theta.Sigma(:,:,j)));
end

% for the missing data points calculate the density over the marginal distribution of the observed components
pointsMissing = find(~pointsObserved);
for n = 1:length(pointsMissing)
    observedComponents = ~isnan(X(pointsMissing(n),:));
    for j = 1:K
        densEst(pointsMissing(n),j) = mvnpdf(X(pointsMissing(n),observedComponents),theta.mu(observedComponents,j)',squeeze(theta.Sigma(observedComponents,observedComponents,j)));
    end
end

end
                                                                                                                                            