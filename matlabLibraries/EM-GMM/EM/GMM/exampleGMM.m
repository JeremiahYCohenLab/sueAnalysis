% set parameters for sampling
mix(1) = 0.3;
mix(2) = 0.7;
mu(:,1) = [3 1];
mu(:,2) = [0 3];
Sigma(:,:,1) = [1 0; 0 1];
Sigma(:,:,2) = [2 1; 1 2];

% sample from the mixture: sample z, then sample x
n = 100;
z = (rand(n,1)>mix(1))+1;
X = zeros(n,2);
for i = 1:2
	X(z==i,:) = mvnrnd(mu(:,i)',squeeze(Sigma(:,:,i)),sum(z==i));
end

% plot true model
GMM_plot(X,mix,mu,Sigma);

% set X values randomly to zero
XNan = X;
ratioNan = 0.3;
XNan(rand(size(X))<ratioNan) = NaN;
% remove rows that have only NaN values
XNan(sum(isnan(XNan),2)==2,:) = [];

% run the EM-algorithm
[expectations theta] = fitGaussMixture(X,2,'kmeans');
GMM_plot(X,theta.tau,theta.mu,theta.Sigma);

[expectationsNan thetaNan] = fitGaussMixture(XNan,2,'kmeans');
GMM_plot(X,thetaNan.tau,thetaNan.mu,thetaNan.Sigma);

