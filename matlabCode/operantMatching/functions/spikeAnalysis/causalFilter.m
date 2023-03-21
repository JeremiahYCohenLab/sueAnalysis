function kernel = causalFilter(len, tauR, tauD, delay)
decay = exp(-(1:len)/tauD);
rise = 1 - exp(-(1:len)/tauR);
kernel = [zeros(1, delay) decay.*rise];