function [p, csF, csFT] = timeProjection(ledLL, csT, startFrame, startTrial, ratio)
    csFT = [];     
    for j = 1:length(startFrame)
        if j < length(startFrame)
        FTTemp = round(ratio/1000*(csT(startTrial(j): startTrial(j + 1)-1) - csT(startTrial(j))) + startFrame(j));
        else
            FTTemp = round(ratio/1000*(csT(startTrial(j): end) - csT(startTrial(j))) + startFrame(j));
        end           
        csFT = [csFT, FTTemp];
    end
    csFT = csFT(csFT < length(ledLL));
    csF = zeros(1,length(ledLL));
    csF(csFT) = 1;
    cskernel = ones(1, round(0.5 * ratio));
    csF = conv(csF, cskernel);
    csF = csF(1:(end - length(cskernel) + 1));
    p = (2*csF) * ledLL;
end
