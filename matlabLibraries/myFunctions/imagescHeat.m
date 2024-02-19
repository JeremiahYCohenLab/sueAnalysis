function h = imagescHeat(image, thresh)
imagesc(image);

maxV = max(image, [], 'all');
minV = min(image, [], 'all');

Nstep = 1000;
Pstep = round(Nstep*(maxV - thresh)/(thresh - minV));
myMap = [[linspace(0, 1, Nstep)' linspace(0, 1, Nstep)' ones(Nstep, 1)];...
[ones(Pstep, 1) linspace(1, 0, Pstep)' linspace(1, 0, Pstep)']];

colormap(gca, myMap)
colorbar
clim([minV maxV])



