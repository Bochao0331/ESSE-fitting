function out = mixed_rho_T2_T2prime_B0_fitting(x,echoTime,y,TEs,zero_shift_ESSE)
% -----------------------
% Input:
%   x(1):   real part of rho
%   x(2):   imaginary part of rho
%   x(3):   R2 (1/msec)
%   x(4):   R2' (1/msec)
%   x(5): off-reasoance (Hz)
%   echoTime:    vector of echo time
%   y:      obeserved data (2xN,1)
%   TE:     vecotr of spin-echo time
%   zero_shift_ESSE: zero-shifted echo time in ESSE

real_rho = x(1);
imag_rho = x(2);
R2 = x(3);
R2prime = x(4);
delta_f = x(5)/1000 ;

Ax = zeros(length(y),1);

for nTE = 1:length(echoTime) % construct vector of signal at echa time from signal equations
    te = echoTime(nTE);
    if any(te == TEs)
        Ax(nTE) = real_rho * exp(-R2*te) ;
        Ax(nTE+length(echoTime)) = imag_rho * exp(-R2*te);
    else
        if  te < zero_shift_ESSE
            Ax(nTE) = ( real_rho * cos(2*pi*delta_f*(te-zero_shift_ESSE)) - imag_rho * sin(2*pi*delta_f*(te-zero_shift_ESSE)) ) * exp(-te*(R2-R2prime)) * exp(-R2prime*zero_shift_ESSE);
            Ax(nTE+length(echoTime)) = ( imag_rho * cos(2*pi*delta_f*(te-zero_shift_ESSE)) + real_rho *sin(2*pi*delta_f*(te-zero_shift_ESSE)) ) * exp(-te*(R2-R2prime)) * exp(-R2prime*zero_shift_ESSE);
        elseif te > zero_shift_ESSE
            Ax(nTE) = ( real_rho * cos(2*pi*delta_f*(te-zero_shift_ESSE)) - imag_rho * sin(2*pi*delta_f*(te-zero_shift_ESSE)) ) * exp(-te*(R2+R2prime)) * exp(R2prime*zero_shift_ESSE);
            Ax(nTE+length(echoTime)) = ( imag_rho * cos(2*pi*delta_f*(te-zero_shift_ESSE)) + real_rho *sin(2*pi*delta_f*(te-zero_shift_ESSE)) ) * exp(-te*(R2+R2prime)) * exp(R2prime*zero_shift_ESSE);
        end
    end
end

out = Ax - y;
end
