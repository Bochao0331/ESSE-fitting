function out = rho_T2_T2prime_B0_fitting(x,TEs,y,TE)
% -----------------------
% Input:
%   x(1):   real part of rho
%   x(2):   imaginary part of rho
%   x(3):   R2 (1/msec)
%   x(4):   R2' (1/msec)
%   x(5): off-reasoance (Hz)
%   TEs:    vector of echo-shift time
%   y:      obeserved data (2xN,1)
%   TE:     spin-echo time

real_rho = x(1);
imag_rho = x(2);
R2 = x(3);
R2prime = x(4);
delta_f = x(5)/1000 ;
  
Ax = zeros(length(y),1);
% for nTE = 1:length(TEs) % construct vector of signal at echa time from signal equations
%     te = TEs(nTE);
%     if  te < TE
%         Ax(nTE) = ( real_rho * cos(2*pi*delta_f*(te-TE)) - imag_rho * sin(2*pi*delta_f*(te-TE)) ) * exp(-(te-TE)*(R2-R2prime));
%         Ax(nTE+length(TEs)) = ( imag_rho * cos(2*pi*delta_f*(te-TE)) + real_rho *sin(2*pi*delta_f*(te-TE)) ) * exp(-(te-TE)*(R2-R2prime));
%     elseif te == TE
%         Ax(nTE) = real_rho ;
%         Ax(nTE+length(TEs)) = imag_rho ;
%     elseif te > TE
%         Ax(nTE) = ( real_rho * cos(2*pi*delta_f*(te-TE)) - imag_rho * sin(2*pi*delta_f*(te-TE)) ) * exp(-(te-TE)*(R2+R2prime));
%         Ax(nTE+length(TEs)) = ( imag_rho * cos(2*pi*delta_f*(te-TE)) + real_rho *sin(2*pi*delta_f*(te-TE)) ) * exp(-(te-TE)*(R2+R2prime));
%     end
% end

for nTE = 1:length(TEs) % construct vector of signal at echa time from signal equations
    te = TEs(nTE);
    if  te < TE
        Ax(nTE) = ( real_rho * cos(2*pi*delta_f*(te-TE)) - imag_rho * sin(2*pi*delta_f*(te-TE)) ) * exp(-te*(R2-R2prime)) * exp(-R2prime*TE);
        Ax(nTE+length(TEs)) = ( imag_rho * cos(2*pi*delta_f*(te-TE)) + real_rho *sin(2*pi*delta_f*(te-TE)) ) * exp(-te*(R2-R2prime)) * exp(-R2prime*TE);
    elseif te == TE
        Ax(nTE) = real_rho * exp(-R2*TE) ;
        Ax(nTE+length(TEs)) = imag_rho * exp(-R2*TE);
    elseif te > TE
        Ax(nTE) = ( real_rho * cos(2*pi*delta_f*(te-TE)) - imag_rho * sin(2*pi*delta_f*(te-TE)) ) * exp(-te*(R2+R2prime)) * exp(R2prime*TE);
        Ax(nTE+length(TEs)) = ( imag_rho * cos(2*pi*delta_f*(te-TE)) + real_rho *sin(2*pi*delta_f*(te-TE)) ) * exp(-te*(R2+R2prime)) * exp(R2prime*TE);
    end
end

out = Ax - y;
end
