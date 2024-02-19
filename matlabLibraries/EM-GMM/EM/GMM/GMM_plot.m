function GMM_plot(X,tau,mu,Sigma)


diffX = max(X(:,1))-min(X(:,1));
diffY = max(X(:,2))-min(X(:,2));
rangeX = [min(X(:,1))-diffX/10 max(X(:,1))+diffX/10];
rangeY = [min(X(:,2))-diffY/10 max(X(:,2))+diffY/10];

vecX = linspace(rangeX(1),rangeX(2),50);
vecY = linspace(rangeY(1),rangeY(2),50);
[XX YY] = meshgrid(vecX,vecY); XX = reshape(XX,1,[]); YY = reshape(YY,1,[]);

figure
subplot(2,1,1)
scatter(X(:,1),X(:,2)); hold on;
for i = 1:length(tau)
    evalF = mvnpdf([XX; YY]',mu(:,i)',Sigma(:,:,i));
    contour(vecX,vecY,reshape(evalF,length(vecX),[]));
end
box on
title('Class-conditionals');


subplot(2,1,2)
evalF = zeros(length(XX),1);
for i = 1:length(tau)
	evalF = evalF + tau(i)*mvnpdf([XX; YY]',mu(:,i)',Sigma(:,:,i));
end

for i = 1:length(tau)
	z(:,i) = tau(i)*mvnpdf(X,mu(:,i)',Sigma(:,:,i));
end
z = z./repmat(sum(z,2),1,2);
[~,idx] = max(z,[],2);	
h = gscatter(X(:,1),X(:,2),idx,[1 0.7 0; 0.85 0.3 0.3],'os',5); hold on;
for i = 1:length(tau)
	set(h(i),'MarkerFaceColor',get(h(i),'Color'));
end
contour(vecX,vecY,reshape(evalF,length(vecX),[]));
set(gca,'XLim',[rangeX(1) rangeX(2)],'YLim',[rangeY(1) rangeY(2)]);
legend off
title('Mixture density');

