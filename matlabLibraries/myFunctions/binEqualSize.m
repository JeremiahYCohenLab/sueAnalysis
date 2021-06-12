function edges = binEqualSize(x, numBins)
% result in numBins or numBins + 1
% edges includes left edge, last one includes right edge
    [~, I] = sort(x);
    binSize = round(length(x)/numBins);
   % extra = floor(0.5 * rem(length(x),numBins));
    edges = min(x);
    edgePosi = 0;
    endFlag = false;
    for b = 1:numBins
        edgePosi = edgePosi+(binSize-1);
        if edgePosi+1 >= length(x)
            endFlag = true;
            edges = [edges, max(x)];
            break
        end
        while x(I(edgePosi)) == x(I(edgePosi+1)) && ~endFlag
            edgePosi = edgePosi + 1;
            if edgePosi + 1 >= length(x)
                endFlag = true;
            end
        end
        if ~endFlag 
            edges = [edges, x(I(edgePosi+1))];
        else
            edges = [edges, max(x)];
            break
        end
    end
    
    if ~endFlag
        if length(x)-edgePosi>0.5*binSize
            edges = [edges, max(x)];
        else
            edges(end) = max(x);
        end
    end
end