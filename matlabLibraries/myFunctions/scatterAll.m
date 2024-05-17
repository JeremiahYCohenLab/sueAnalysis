function scatterAll(matrix, names, spotSize, color)
% matrix: n * x matrix, n:sample number, x:features to scatter
% names: names of the features, 1* n cell of strings

s = size(matrix,2);
for i = 1:s
    for j = i:s
        if i == j
            subplot(s, s,(s)*(i-1)+j); hold on;
            histogram(matrix(:,i), 50, 'FaceColor', color, 'EdgeColor','none')
            title(names{i})
        else
             subplot(s,s,(s)*(i-1)+j); hold on;
             scatter(matrix(:,i), matrix(:,j), spotSize, color, 'filled');
             line(minmax(matrix(:,i)'), [0 0],  'color', [0.7 0.7 0.7], 'LineStyle','--')
             line([0 0], minmax(matrix(:,j)'), 'color', [0.7 0.7 0.7], 'LineStyle','--')    
             xlabel(names{i})
             ylabel(names{j})
             valInd = ~isnan(matrix(:,i))& ~isnan(matrix(:,j));
             [R, P] = corrcoef(matrix(valInd,i), matrix(valInd,j));
             title(sprintf('%0.2f p:%0.2f', R(1,2), P(1,2)));
        end
    end
end