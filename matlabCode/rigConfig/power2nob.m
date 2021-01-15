function nob = power2nob(power)
% map laser to powernob number
% rig 295F neuralynx
% last recalibration: 08/18/2020 by ZS
 nob = log((power + 14.83)/4.365)/0.3459; 