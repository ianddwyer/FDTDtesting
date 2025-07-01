function [kappa, bcoef, ccoef] = ...
    PmlScaleEUpdate(iend,del,dt,eps,npml,m,ma,sigmax,kapmax,alpmax)
% Function used to compute the PML coefficient arrays used for
% e-field updates in the PML region.
%------------------------------------------------------------
% Input:
% iend = 1 (left side), or = 2 (right side)
% npml = the number of pml cells
% m = polynomial scaling for sigma and kappa
% ma = polynomial scaling for alpha
% sigmax,kapmax,alpmax = max values for sigma, kappa, and alpha
%---------------------------------------------------------------
sigma = zeros(npml+1,1);
kappa = ones(npml+1,1);
alpha = zeros(npml+1,1);
bcoef = zeros(npml+1,1);
ccoef = zeros(npml+1,1);
if(npml < 1) 
    return;
end
if(iend == 1)
    for i = 1:npml+1
        polyscalesk = (double(npml+1-i)/double(npml))^m;
        polyscalea = (double(i -1)/double(npml))^ma;
        sigma(i) = sigmax * polyscalesk;
        kappa(i) = 1.0 + double(kapmax-1.0) * polyscalesk;
        alpha(i) = alpmax * polyscalea;
    end
else
    for i = 1:npml+1
        polyscalesk = (double(i -1)/double(npml))^m;
        polyscalea = (double(npml+1-i)/double(npml))^ma;
        sigma(i) = sigmax* polyscalesk;
        kappa(i) = 1.0 + double(kapmax-1.0)* polyscalesk;
        alpha(i) = alpmax * polyscalea;
    end
end

for i = 1:npml+1
    temp = kappa(i)*alpha(i)+sigma(i);
    dtotau = temp*dt/(kappa(i)*eps);
    bcoef(i) = exp(-dtotau);
    if(temp == 0.0)
        ccoef(i) = 0.0;
    else
        ccoef(i) = sigma(i)*(bcoef(i)-1.0)/(kappa(i)*temp*del);
    end
end