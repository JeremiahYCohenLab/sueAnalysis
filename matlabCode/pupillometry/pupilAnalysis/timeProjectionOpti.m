function [ratioMax, p] = timeProjectionOpti(ledLL, csT, startFrame, startTrial, ratio)
    parfor i = 1:length(ratio)
        csFT = [];     
        for j = 1:length(startFrame)
            if j < length(startFrame)
            FTTemp = round(ratio(i)/1000*(csT(startTrial(j): startTrial(j + 1)-1) - csT(startTrial(j))) + startFrame(j));
            else
                FTTemp = round(ratio(i)/1000*(csT(startTrial(j): end) - csT(startTrial(j))) + startFrame(j));
            end           
            csFT = [csFT, FTTemp];
        end
        
        csFT = csFT(csFT < length(ledLL));
        csF = zeros(1,length(ledLL));
        csF(csFT) = 1;
        cskernel = ones(1, round(0.5 * ratio(i)));
        csF = conv(csF, cskernel);
        csF = csF(1:(end - length(cskernel) + 1));
        p(i) = csF * ledLL;
    end
    [~, m] = max(p);
    ratioMax = ratio(m);
end
