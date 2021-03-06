function [xJumpVec,yJumpVec,totalMagJumpVec] = cFactor(g,AH,params,phyGrid)
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

k = (-floor((N1-1)/2):ceil((N1-1)/2)).';
Sipi = 1.85193705198247;          % normalizing constant
sig = pi*sin( pi*abs(k)/max(k) )/Sipi;

% Test Stuff
% fftShiftMat = zeros(N1);
% for ii = 1:N1
%     basisVec = zeros(N1,1);
%     basisVec(ii) = 1;
%     fftShiftMat(:,ii) = fftshift(basisVec);
% end
% 
% xJumpTest = fftShiftMat*diag(1i*sign(k).*sig)*1/N1*fftShiftMat*mmv_t;


nPts = N1;
h = 2*pi/nPts;
phyGrid = -pi + h*(0:nPts-1).';

if strcmp(params.signal,'real')
    if N2 <= 1
        g_hat_x = 1/N1*fftshift(mmv_t);
    else
        g_hat_x = ifft(1/(N1)*fftshift(mmv_t,1),[],2);
        g_hat_y = ifft(1/(N2)*fftshift(mmv_t.',1),[],2);
    end
elseif strcmp(params.signal,'complex')
    if N2 < 2
        g_hat = 1/N1*conv(mmv_t,repmat(conj(flipud(mmv_t)),3,1),'same');
        g_hat_x = 1/N1*g_hat;
%         g_hat = 1/N1*conv(mmv_t,repmat(conj(flipud(mmv_t)),3,1),'same');
    else
        g_hat = rot90(fftshift(1/(N1*N2)*conv2(mmv_t,repmat(conj(rot90(mmv_t,2)),3,3),'same')),2);
        g_hat_x = ifft(1/(N1)*fftshift(g_hat,1),[],2);
        g_hat_y = ifft(1/(N2)*fftshift(g_hat.',1),[],2);
%         g_hat = 1/(N1*N2)*fftshift(conv2(mmv_t,repmat(conj(rot90(mmv_t,2)),3,3),'same'));
    end
end

% Create trig concentration factor
jump_x_mat = zeros(N1,N2);
fourKern = exp( 1i*phyGrid*k.' );

% g_hat_x = ifft(g_hat_x,[],2);
for jj = 1:N2
    sigTrig = g_hat_x(:,jj).*(1i*sign(k).*sig); 
    jump_x_mat(:,jj) = fftshift(real( fourKern*sigTrig ));
    jump_x_matim(:,jj) = fftshift(imag( fourKern*sigTrig ));
end

% jump_x_mat = real(reshape(AH(fftshift(c_fact.*g_hat)),N1,N2));
xJumpVec = reshape(jump_x_mat,N1*N2,1);

yJumpVec = zeros(N1*N2,1);
if N2 > 1
    jump_y_mat = zeros(N1,N2);
    for jj = 1:N1
        sigTrig = g_hat_y(:,jj).*(1i*sign(k).*sig); 
        jump_y_mat(:,jj) = fftshift(real( fourKern*sigTrig ));
    end
    yJumpVec = reshape(jump_y_mat.',N1*N2,1);
end

if strcmp(params.signal,'complex')
    xJumpVec = xJumpVec/2;
    yJumpVec = yJumpVec/2;
end

totalMagJumpVec = abs(xJumpVec) + abs(yJumpVec);
end