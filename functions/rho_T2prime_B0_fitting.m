function out = rho_T2prime_B0_fitting(x, TEs, y, TE, T2)
% -----------------------
% Input:
%   x(1):   real part of rho
%   x(2):   imaginary part of rho
%   x(3):   T2' (msec)
%   x(4): off-reasoance (Hz)
%   TEs:    vector of echo-shift time
%   y:      obeserved data (2xN,1)
%   TE:     spin-echo time

real_rho = x(1);
imag_rho = x(2);
R2prime = x(3);
delta_f = x(4)/1000 ;
  
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
