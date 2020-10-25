function laserInput = desiredLaserInput_295acute(desiredOutput, laserSource, opticFiberTrans)
% laserInput = desiredLaserInput(desiredOutput, opticFiberTrans)
%   DESCRIPTION
%       Uses a calibration range from 1 to 10 fit with output = a*exp(b*input)
%   OUTPUT(S)
%       laserInput: knob setting for a given laser to achive a desired output
%   INPUT(S)
%       desiredOutput: output at most distal optic fiber (mW)
%       laserSource: integer; either '532' or '473'
%       opticFiberTrans: vector of transmission percentage for optic
%       fiber(s)96
if nargin < 3
    opticFiberTrans = 0.8;
end

% both calibrated on 20181010 for split patch cords
% 473 calibrated on 20190930 (still split patch cord)
if laserSource == 473
   laserInput = log(85.28/(desiredOutput*2/opticFiberTrans) - 1) / -0.9716 + 7.902;
elseif laserSource == 532
    p = [0.295 10.73 (-23.84 - desiredOutput*2/opticFiberTrans)];
    r = roots(p);
    laserInput = max(r);
else
    error('laserSource must be either 473 or 532')
end