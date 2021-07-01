function [img_List,c_img, kdata, info] = readESSE_twix(file, flags)
%   Read multiple TSE dataset from TWIX raw data
%   Input:
%       - file: path of raw data files

%   Output:
%   Saved as cell format, each cell contains one echo-shift k-space or
%   images
%       - img_List: cell of coil-combined images (Nx x Ny x Nslice)
%       - c_img: cell of coil images (Nx x Ny x Ncoil x Nslice)
%       - kdata: cell of coil k-space (Nx x Ny x NCoil x Nslice)
%       - info: basic informatino of raw data

%   Author: Bochao Li
%   Email: bochaoli@usc.edu

twix_Files = dir(file);
nF = length(twix_Files);

for nF = 1:length(twix_Files)
    baseFileName = twix_Files(nF).name;
    fullFileName = fullfile(flags.folder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    twix = mapVBVD(fullFileName);
    
    if length(twix) ==2
        twix = twix{2};
    end
    info.TurboFactor = twix.image.NSeg;
    
    % ----Set twix flag----
    twix.image.flagIgnoreSeg = flags.IgnoreSeg;
    twix.image.flagDoAverage = flags.DoAverage;
    os  = flags.os; % oversampling factor
    
    % ----Reconstruction----
    img = zeros([twix.image.NCol/os, twix.image.NLin, twix.image.NSli], 'single'); %sum-of-square coil-combined images
    for sli = 1:twix.image.NSli
        %--read in the data for slice 'sli'--
        temp_kdata = permute(twix.image(:,:,:,1,sli),[1 3 2]);
        
        %----------------------------------------------------------------------
        % Remove oversampling
        %----------------------------------------------------------------------
        if flags.RemoveOS == true
            Nkx = size(temp_kdata,1);
            Nx = Nkx/2;
            imc_ = fftshift(fft(ifftshift(temp_kdata, 1), [], 1), 1); % Nkx x Nky x Ncoil
            idx1_range = (-floor(Nx/2):ceil(Nx/2)-1).' + floor(Nkx/2) + 1;
            imc_ = imc_(idx1_range,:,:,:);  % Nx x Nky x Nkz x Nc
            temp_kdata = fftshift(ifft(ifftshift(imc_, 1), [], 1), 1); % Nkx x Nky x Ncoil           
        end
        kdata{nF}(:,:,:,sli) = temp_kdata;
        c_img{nF}(:,:,:,sli) = fft2c(temp_kdata); % coil images
        img(:,:,sli) = squeeze(sqrt(sum(abs(c_img{nF}(:,:,:,sli)).^2,3))); % sum-of-square coil-combined images (maybe abandon later)
    end
    img_List{nF} = img; % save each echo-shift images
end
    info.TE = twix.hdr.MeasYaps.alTE{1}; % TE: usec
    info.TR = twix.hdr.MeasYaps.alTR{1}; % TR: usec
end

