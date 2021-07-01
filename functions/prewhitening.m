function [kdata_prewhitening] = prewhitening(kdata,inv_L)
%   Calculae noise covariance marix from noise samples
%   Input:
%       - kdata: k-space data (Nx x Ny x Nslice x Ncoil)

%   Output:
%       - Psi: scaled noise coavariance matrix (Ncoil x Ncoil)
%       - inv_L: inverse sqaure roof of the noise covariance (Ncoil x Ncoil)

%   Author: Bochao Li
%   Email: bochaoli@usc.edu

[Nx, Ny, Nslice, Ncoil] = size(kdata);
kdata_prewhitening = ipermute(reshape(inv_L * reshape(permute(kdata,[4 1 2 3]), [Ncoil Nx*Ny*Nslice]), [Ncoil Nx Ny Nslice]), [4 1 2 3]);

end

