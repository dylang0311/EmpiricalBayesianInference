function output = calc_ci(output,BI,signal,sigT)
% calculate 95% confidence interval of output.x given burn in BI, signal
% type, and transform to signal domain sigT

[N,~] = size(output.x); 

output.ci = zeros(2,N);
output.x = sigT(output.x);

for ii = 1:N 
    if strcmp(signal,'real')
        dat = output.x(ii,BI:end); 
    elseif strcmp(signal,'complex')
        dat = abs(output.x(ii,BI:end)); 
        
        datAng = angle(output.x(ii,BI:end)); 
        output.ciAng(1,ii) = quantile(datAng,0.025);
        output.ciAng(2,ii) = quantile(datAng,0.975);
    end
    output.ci(1,ii) = quantile(dat,0.025);
    output.ci(2,ii) = quantile(dat,0.975);
end


end

