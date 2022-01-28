function [xJumpVec,yJumpVec,totalMagJumpVec] = cFactor(g,AH,params)
% Uses the method of concentration factors to find the jump function
% directly from Fourier data.
% 
% Inputs:
% g = input signal
% AH = Inverse of forward operator
% params.signal = 'real' or 'compex'
% params.N1,params.N2 = dimensions
%
% Outputs:
% xJumpVec = jump function in first dimension
% yJumpVec = jump function in second dimension
% totalMagJumpVec = absolute value of total jump function

N1 = params.N1;
N2 = params.N2;

mmv_t = reshape(g,N1,N2);
if strcmp(params.signal,'real')
    g_hat = fftshift(mmv_t);
elseif strcmp(params.signal,'complex')
    if N2 < 2
        g_hat = 1/N1*conv(mmv_t,repmat(conj(flipud(mmv_t)),3,1),'same');
    else
        g_hat = 1/(N1*N2)*fftshift(conv2(mmv_t,repmat(conj(rot90(mmv_t,2)),3,3),'same'));
    end
end

% Create trig concentration factor
k = (-(N1-1)/2:(N1-1)/2)';
Sipi = 1.85193705198247;          % normalizing constant
sig = pi*sin( pi*abs(k)/max(k) )/Sipi;
c_fact = (1i*sign(k).*sig);

jump_x_mat = abs(reshape(AH(c_fact.*g_hat),N1,N2));
xJumpVec = reshape(jump_x_mat,N1*N2,1);

yJumpVec = zeros(N1*N2,1);
if N2 > 1
    jump_y_mat = abs(reshape(AH(g_hat.*c_fact'),N1,N2));
    yJumpVec = reshape(jump_y_mat,N1*N2,1);
end

if strcmp(params.signal,'complex')
    xJumpVec = xJumpVec/2;
    yJumpVec = yJumpVec/2;
end

totalMagJumpVec = abs(xJumpVec) + abs(yJumpVec);
end