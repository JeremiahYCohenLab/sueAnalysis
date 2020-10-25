function [startFrame, ratio] = timeOpti(session, iter, errorMax)
errorSize = 1;
errorNum = 100;
for i = 1:20
    [cs, startFrame, ratioMax, ] = timeAlign(session, plotFlag);
    perf = round(ratioMax * 0.5);
    errorNum = length(find(cs < perf - 2*errorSize));
    if errorNum <= errorMax
        print('Criterion reached')
        break;
    end
    
end

if i == iter
    print('Maximum Iter reached')
end
end
