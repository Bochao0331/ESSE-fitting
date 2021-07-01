%   This scrip jointly estimates pronton density, T2,T2',off-resonance 
%   from 2D echo-sfhited SE data with different echo-sfhited time (non-linear least-sqaure fitting)

%   1. You can read either ISMRMD or MAT(generated based on ISMRMD) fiels.

%	2. For TWIX Raw data, the order of acquired data: negative echo-shift
%  --> zero-shift --> postive-echo shift (Check the ISMRMD files' names for details)

%   3. MAT file contains cells (each cells contains one echo-shift data; same echo-shift order as TWIX):
%       (1) k-space date: k_data_all (Nx x Ny x Ncoil x Nz)
%       (2) coil-combined images (magnitude,sum-of-square combination): img_List_all (Nx x Ny x Nslice)
%       (3) coil images: c_img_all (Nx x Ny x Ncoil x Nslice)Nslice
%       (4) info: information obtained from TWIX (TE,TR,TurboFactor...)

%   Example Raw data locate in: /mnt/sdata_new/Bochao/0528AppleESSE/RAW/ on MREL
%   server; 
%   Example MAT files locate in: /mnt/sdata_new/Bochao/0528AppleESSE//MAT

%	Author: Bochao Li
%	Email: bochaoli@usc.edu

clear
%% add path
%% add path
addpath(genpath('.\functions'));
computer_type = computer;
if strcmp(computer_type, 'PCWIN64')
    ismrmrd_directory = 'F:\USC\MREL\Tool\ISMRMD';
elseif strcmp(computer_type, 'GLNXA64')
    src_directory = '';
    ismrmrd_directory = '/server/home/nlee/ismrmrd';
    data_parent_directory = '';
end
addpath(genpath(ismrmrd_directory));

%% ---- Reading parameter settings ----
zero_shift_ind = 5; % zero-shift data index
data_format = 'ismrmd'; % Data reading option 'twix'(Raw) or 'mat';
Mat_folder = 'F:\USC\MREL\LowField\LungImaging\T2measurement\Data\0528AppleGroundTruth\MAT'; % Direcoty where MAT saved
Mat_file = fullfile(Mat_folder,'ES_SE.mat'); % MAT name

%% ---- Loading data ----
if strcmp(data_format, 'ismrmd')
    h5_folder = 'F:\USC\MREL\LowField\LungImaging\T2measurement\Data\0528AppleGroundTruth\RAW\h5';
    noise_folder = 'F:\USC\MREL\LowField\LungImaging\T2measurement\Data\0528AppleGroundTruth\RAW\noise';
    h5_fileList = fullfile(h5_folder, '*tse_ES_BLi*.h5'); % Change to whatever pattern you need.
    noise_fileList = fullfile(noise_folder, '*tse_ES_BLi*.h5'); % Change to whatever pattern you need.
    h5_Files = dir(h5_fileList);
    noise_Files = dir(noise_fileList);
    nFiles = length(h5_Files);
    
    
    % ----ISMRMD .h5 reading flags----
    Read_flags.h5_fileList       = h5_fileList;
    Read_flags.noise_fileList    = noise_fileList;
    Read_flags.RemoveOS          = true; % remove oversampling
    Read_flags.IgnoreSeg         = true; % concatanate segmemtns
    Read_flags.DoAverage         = true; % Do averages (if 'average' was set during data acquistion)
    Read_flags.CropPhaseEncoding = true;
    Read_flags.Squeeze           = true;
    Read_flags.os                = 2; % oversampling rate (Siemens Default value, don't change it)
    Read_flags.noisePreWhitening = true;
    
    % ----Loading Raw Data----
     [c_img_all, kdata_all, noise_all, info_all] = readSE_ismrmd(Read_flags); %loading data from TWIX.
    if ~isfile(Mat_file)
        save (Mat_file, 'c_img_all', 'kdata_all','noise_all','info_all');
    end
    
elseif strcmp(data_format, 'mat')
    load(Mat_file);
end

%% ---- Get data info ----
[Nx, Ny, Nslice, Ncoil] = size(kdata_all{1});
TE = zeros(length(info_all), 1,'double');
for nF = 1:length(info_all)
    TE(nF) = info_all{nF}.TE; % msec
end
echoshifts = [-12 -5 -2.3  -1 0 1 2.3 5 12].'; % echo shift value
echo_time_all = TE + echoshifts; % Echo time of the echo-shift

%% Apply the noise prewhitening matrix on k-space before Recon
% ---- Calculate noise covariance ----
for k = 1:length(noise_all)
    [Psi, inv_L] = calculate_noise_covariance(noise_all{k});
    kdata_prew{k} = prewhitening(kdata_all{k},inv_L);
end

%% ---- Coil Sensitivity Estimation & Coil-comined Recon ----
csm_option.method = 'walsh'; % coil estimation method: 'walsh', 'sos'
csm_option.cal_shape = [32 32]; %calibration region
csm_option.kdata = kdata_prew{zero_shift_ind};
[csm, cal_im] = coil_estimation(csm_option);

load("csm_se.mat");
% ---- Coil combination ----
for nD = 1 : length(c_img_all)
    cc_img_List{nD} = squeeze(sum( conj(csm) .* c_img_all{nD}, 4)); % Coil combined images (Nx x Ny x Nslice for each echo-shift cell)
end


%% Show Recon images
fSIZE = 16;
XTICK = [ size(cc_img_List{1},2) * [0.5:1:size(cc_img_List,2)] ];

[ha, pos] = tight_subplot(4,1,[.01 .01],[.04 .01],[.02 .03]);
axes(ha(1));  imagesc(1e2* abs(cell2mat(cc_img_List)), [0 1]); colormap(gca,gray); axis image; colorbar; ylabel('Magnitude')
axes(ha(2));  imagesc(angle(cell2mat(cc_img_List))/pi*180,[-180 180]); colormap(gca,hsv); axis image; cb = colorbar; ylabel('Phase');ylabel(cb,'degree','FontWeight','BOLD')
axes(ha(3));  imagesc(1e2* real(cell2mat(cc_img_List)), [-1 1]); colormap(gca,gray); axis image; colorbar; ylabel('Real')
axes(ha(4));  imagesc(1e2* imag(cell2mat(cc_img_List)), [-1 1]); colormap(gca,gray); axis image; colorbar; ylabel('Imaginary')

set(ha,'XTickLabel','','FontSize',fSIZE,'FontWeight','BOLD'); set(ha,'YTickLabel','')
set(ha(4),'XTick',XTICK, 'XTickLabel',echo_time_all)
set(get(ha(4),'Xlabel'),'String','Echo-shited Time (ms)')

%% ----Generate the ROI mask----
temp = abs(cc_img_List{zero_shift_ind})./mean(mean(abs(cc_img_List{zero_shift_ind}))); % use fist TE images
temp = temp<1; % threshold value for ROI

se = strel('line',4,4);
bw_img = ~imdilate(temp,se); % erosion
bd_bw = bwboundaries(bw_img); % boundary of ROI

% ---- Show ROI ----
figure,
imagesc(abs(cc_img_List{1}));colormap(gray)
hold on;
for k = 1:length(bd_bw)
    boundary = bd_bw{k};
    plot(boundary(:,2), boundary(:,1), 'r--', 'LineWidth', 0.5)
end
axis off image

%% Non-linear fitting (unkonwn T2)
temp = cc_img_List;
for k = 1: length(temp)
    signal_TE(:,k) = 1e2*temp{k}(bw_img); % each columne for one TE; each row for one voxel in ROI
end

tic
parfor k = 1:size(signal_TE,1)
    k
    y = [real(signal_TE(k,:)).'; imag(signal_TE(k,:)).'];
    x0 = [real(signal_TE(k,zero_shift_ind)), imag(signal_TE(k,zero_shift_ind)), 1/50, 1/13, 50]; % inital values
    E = @(x) rho_T2_T2prime_B0_fitting(x,echo_time_all,y, TE(1));
    options = optimoptions('lsqnonlin');
    [x(k,:),resnorm,residual,exitflag,output] = lsqnonlin(E,x0,[],[],options);
end
toc

%% Non-linear fitting (known T2)
load('F:\USC\MREL\LowField\LungImaging\T2measurement\Data\0528AppleGroundTruth\MAT\rho_R2_SE.mat'); % load T2 estimated from SE
T2_array = 1./x(:,3); 
clear x

temp = cc_img_List;
for k = 1: length(temp)
    signal_TE(:,k) = 1e2*temp{k}(bw_img); % each columne for one TE; each row for one voxel in ROI
end

tic
parfor k = 1:size(signal_TE,1)
    k
    T2 = T2_array(k);
    y = [real(signal_TE(k,:)).'; imag(signal_TE(k,:)).'];
    x0 = [real(signal_TE(k,zero_shift_ind)), imag(signal_TE(k,zero_shift_ind)), 1/13, 50]; % inital values
    E = @(x) rho_T2prime_B0_fitting(x,echo_time_all,y, TE(1), T2_array(k));
    options = optimoptions('lsqnonlin');
    [x(k,:),resnorm,residual,exitflag,output] = lsqnonlin(E,x0,[],[],options);
end
toc

%% Non-linear fitting (known  rho & T2)
load('F:\USC\MREL\LowField\LungImaging\T2measurement\Data\0528AppleGroundTruth\MAT\rho_R2_SE.mat'); % load T2 estimated from SE
real_rho_array = x(:,1);
img_rho_array = x(:,2);
T2_array = 1./x(:,3); 
clear x

temp = cc_img_List;
for k = 1: length(temp)
    signal_TE(:,k) = 1e2*temp{k}(bw_img); % each columne for one TE; each row for one voxel in ROI
end

tic
parfor k = 1:size(signal_TE,1)
    k
    T2 = T2_array(k);
    real_rho = real_rho_array(k);
    img_rho = img_rho_array(k);
    y = [real(signal_TE(k,:)).'; imag(signal_TE(k,:)).'];
    x0 = [1/13, 50]; % inital values
    E = @(x) T2prime_B0_fitting(x,echo_time_all,y, TE(1), T2, real_rho, img_rho);
    options = optimoptions('lsqnonlin');
    [x(k,:),resnorm,residual,exitflag,output] = lsqnonlin(E,x0,[],[],options);
end
toc

%% T2 map 
T2map_ESSE = zeros(size(bw_img));
T2map_ESSE(bw_img) = 1./x(:,3);

figure,
imagesc(T2map_ESSE,[0 200])
h = colorbar;
ylabel(h, 'T2 (ms)','FontWeight','bold')
colormap(hot)
axis image off
title('T2 map from ES-SE')

% Mat_file = fullfile(Mat_folder,'T2map_ESSE.mat'); % MAT name
% if ~isfile(Mat_file)
%     save (Mat_file,'T2map_ESSE');
% end

%% T2' map 
T2primemap_ESSE = zeros(size(bw_img));
T2primemap_ESSE(bw_img) = 1./x(:,4);

figure,
imagesc(T2primemap_ESSE,[])
h = colorbar;
ylabel(h, '%','FontWeight','bold')
colormap(jet)
axis image off
% title('T2'' map from ES-SE')

% Mat_file = fullfile(Mat_folder,'T2primemap_ESSE.mat'); % MAT name
% if ~isfile(Mat_file)
%     save (Mat_file,'T2primemap_ESSE');
% end

%% B0 map 
B0map_ESSE = zeros(size(bw_img));
B0map_ESSE(bw_img) = B0_knownT2;

figure,
imagesc(B0map_ESSE,[-60 60])
h = colorbar;
ylabel(h, '\Deltaf (Hz)','FontWeight','bold')
colormap(hsv)
axis image off
title('\Deltaf map from ES-SE')

% Mat_file = fullfile(Mat_folder,'B0map_ESSE.mat'); % MAT name
% if ~isfile(Mat_file)
%     save (Mat_file,'B0map_ESSE');
% end

%% rho map 
rhomap_ESSE = zeros(size(bw_img));
rhomap_ESSE(bw_img) = rho_R2prime_B0_ESSE_TE12csm;

figure,
imagesc(real(rhomap_ESSE),[0 3])
h = colorbar;
% ylabel(h,'%','FontWeight','bold')
colormap(gray)
axis image off
title('\rho real map from rho_R2prime_B0_ESSE_TE12csm')


%% Display a signal or phase from an individual voxel
im_tse = reshape(cell2mat(cc_img_List),[Nx Ny Nslice length(cc_img_List)]);
figure;
subplot(1,2,1); imagesc(abs(im_tse(:,:,:,1))); axis image off
colormap(gray(256));
for idx = 1:1000
    [x,y,button] = ginput(1);
    row = floor(y)
    col = floor(x)
    subplot(1,2,2);
    plot(echo_time_all, reshape(angle(im_tse(row,col,:,:))/pi*180, [length(echo_time_all) 1]),'o-');
    xline(echo_time_all(zero_shift_ind),'r--')
    yline(angle(im_tse(row,col,:,zero_shift_ind))/pi*180,'r--')
    ylabel('Phase(degree)','FontWeight','BOLD')
    xlabel('Echo-shifted Time (ms)','FontWeight','BOLD')
    title('Phase of ES-SE data','FontWeight','BOLD')
    set(gca,'XTick',echo_time_all, 'XTickLabel',echo_time_all)
    if button == 3
        return;
    end
end

%% Display a signal or phase from an individual voxel
im_tse = reshape(cell2mat(cc_img_List),[Nx Ny Nslice length(cc_img_List)]);
figure;
subplot(1,2,1); imagesc(abs(im_tse(:,:,:,1))); axis image off
colormap(gray(256));
for idx = 1:1000
    [x,y,button] = ginput(1);
    row = floor(y)
    col = floor(x)
    subplot(1,2,2);
    plot(echo_time_all, reshape(abs(im_tse(row,col,:,:)), [length(echo_time_all) 1]),'o-');
    xline(echo_time_all(zero_shift_ind),'r--')
    yline(abs(im_tse(row,col,:,zero_shift_ind)),'r--')
    ylabel('Phase(degree)','FontWeight','BOLD')
    xlabel('Echo-shifted Time (ms)','FontWeight','BOLD')
    title('Phase of ES-SE data','FontWeight','BOLD')
    set(gca,'XTick',echo_time_all, 'XTickLabel',echo_time_all)
    if button == 3
        return;
    end
end

%% ----post-analysis T2(T2*) map----

% ---Mean & std in ROI---%
a = mean( T2Upmap(bw_img) ) *1000;
b = std( T2Starmap(bw_img) ) *1000;
%   ----Histogram----%

figure,
for nfk = 1:length(f_k)
    h = histogram(T2map_w{nfk}(bw_img)*1000,10000);
    hold on;
    xlabel('T2*','FontWeight','bold');
    xlim([0 50])
end



%% T2 perncentage difference map 
PerDiffeT2prime = zeros(size(bw_img));
PerDiffeT2prime(bw_img) = PercDiff_T2prime;

figure,
imagesc(PerDiffeT2prime,[-3 3])
h = colorbar;
ylabel(h, 'Percentage difference (%)','FontWeight','bold')
colormap(hot)
axis image off
title('T2'' Percentage difference map from ES-SE')

% Mat_file = fullfile(Mat_folder,'T2map_ESSE.mat'); % MAT name
% if ~isfile(Mat_file)
%     save (Mat_file,'T2map_ESSE');
% end

%% T2 perncentage difference map 
PerDiffeB0 = zeros(size(bw_img));
PerDiffeB0(bw_img) = PercDiff_B0;

figure,
imagesc(PerDiffeB0,[-1 1])
h = colorbar;
ylabel(h, 'Percentage difference (%)','FontWeight','bold')
colormap(hot)
axis image off
title('\Deltaf Percentage difference map from ES-SE')
%% ----exponential fitting (currently we fit seperatly for R2* and R2up(R2-R2'))----
% temp = cc_img_List;
% for k = 1: length(temp)
%     signal_TE(:,k) = temp{k}(bw_img)*1e8; % each columne for one TE; each row for one voxel in ROI
% end
% for k = 1:length(signal_TE)
%     k
%     f = fit(echotime.'/1000,real(signal_TE(k,1:end)).','exp1'); % fitting A*exp(bt)
%     B(k) = -f.b;
%     A(k) = f.a;
% end
% clear signal_TE