function scatterAll(matrix, names, spotSize)
% matrix: n * x matrix, n?sample number, x:features to scatter
% names: names of the features, 1* n cell of strings

s = size(matrix,2);

figure2;
for i = 1:s-1
    for j = i+1:s
     subplot(s-1,s-1,(s-1)*(i-1)+j-1); hold on;
     scatter(matrix(:,i), matrix(:,j), spotSize, 'filled');
     line(minmax(matrix(:,i)'), [0 0],  'color', [0.7 0.7 0.7], 'LineStyle','--')
     line([0 0], minmax(matrix(:,j)'), 'color', [0.7 0.7 0.7], 'LineStyle','--')    
     xlabel(names{i})
     ylabel(names{j})
     [R, P] = corrcoef(matrix(:,i), matrix(:,j));
     title(sprintf('%0.2f p:%0.2f', R(1,2), P(1,2)));
    end
end