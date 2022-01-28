function jumps = freqTV(g,AH,params)
% Inputs:
% g = input signal
% AH = Inverse of forward operator
% params.signal = 'real' or 'compex'
% params.N1,params.N2 = dimensions
%
% Outputs:
% jumps = absolute value of total jump function

[~,~,jumps] = cFactor(g,AH,params);
end