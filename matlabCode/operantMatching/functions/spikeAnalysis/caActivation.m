function signal = caActivation(input, cHalf, slope, peak)
    signal = peak ./(1.+exp(-slope * (input - cHalf)));
end