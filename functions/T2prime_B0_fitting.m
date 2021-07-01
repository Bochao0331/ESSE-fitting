function out = T2prime_B0_fitting(x, TEs, y, TE, T2, real_rho, imag_rho)
% Estimate R2' and off-resonance with given rho and T2 from SE-SE
% -----------------------
% Input:
%   x(1):   R2' (msec)
%   x(2):   off-reasoance (Hz)
%   TEs:    vector of echo-shift time
%   y:      obeserved data (2xN,1)
%   TE:     spin-echo time

R2prime = x(1);
delta_f = x(2)/1000 ;
  
Ax = zeros(length(y),1);
for nTE = 1:length(TEs) % construct vector of signal at echa time from signal equations
    te = TEs(nTE);
    if  te < TE
        Ax(nTE) = ( real_rho * cos(2*pi*delta_f*(te-TE)) - imag_rho * sin(2*pi*delta_f*(te-TE)) ) * exp(-(te-TE)*(1/T2-R2prime));
        Ax(nTE+length(TEs)) = ( imag_rho * cos(2*pi*delta_f*(te-TE)) + real_rho *sin(2*pi*delta_f*(te-TE)) ) * exp(-(te-TE)*(1/T2-R2prime));
    elseif te == TE
        Ax(nTE) = real_rho ;
        Ax(nTE+length(TEs)) = imag_rho ;
    elseif te > TE
        Ax(nTE) = ( real_rho * cos(2*pi*delta_f*(te-TE)) - imag_rho * sin(2*pi*delta_f*(te-TE)) ) * exp(-(te-TE)*(1/T2+R2prime));
        Ax(nTE+length(TEs)) = ( imag_rho * cos(2*pi*delta_f*(te-TE)) + real_rho *sin(2*pi*delta_f*(te-TE)) ) * exp(-(te-TE)*(1/T2+R2prime));
    end
end

out = Ax - y;
end
