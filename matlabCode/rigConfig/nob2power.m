function power = nob2power(nob)
% map nob number to laser power
% rig 295F neuralynx
% last recalibration: 08/18/2020 by ZS
power = 4.365*exp(0.3459*nob) - 14.83;
