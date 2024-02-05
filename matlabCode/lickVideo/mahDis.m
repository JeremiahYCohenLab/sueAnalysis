function dis = mahDis(mat, varargin)
    p = inputParser;
    p.addParameter('cov', []);
    p.parse(varargin{:})

    % calculate non-nans
    ind = ~isnan(sum(mat, 2));
    mat = mat(ind,:);
    m = mean(mat, 1);
    % calculate cov
    if isempty(p.Results.cov) 
        C = cov(mat);
    else
        C = p.Results.cov;
    end
    % calculate distance
    dis = sqrt(diag((mat - m)/(C) * (mat - m)')); 
end