%   This scrip jointly estimates pronton density, T2,T2',off-resonance  from
%   2D echo-sfhited SE data with different echo-sfhited time, and multiple SE data with different TE (non-linear least-sqaure fitting)

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

%% ---- add path ----
addpath(genpath('.\functions'));
SE_path = 'F:\USC\MREL\LowField\LungImaging\T2measurement\Data\0528AppleGroundTruth\MAT\SE_ccImg.mat';
ESSE_path = 'F:\USC\MREL\LowField\LungImaging\T2measurement\Data\0528AppleGroundTruth\MAT\ESSE_ccImg.mat';
SE_data = load(SE_path);
ESSE_data = load(ESSE_path);

SE_ccImg = SE_data.cc_img_List;
SE_info = SE_data.info_all;
ESSE_ccImg = ESSE_data.cc_img_List;
ESSE_info = ESSE_data.info_all;
cc_img_List = {SE_ccImg{1:end}, ESSE_ccImg{1:end}};

%% ---- Get data info of SE----
[Nx, Ny, Nslice] = size(SE_ccImg{1});
SE_TE = zeros(length(SE_info), 1,'double');
for nF = 1:length(SE_info)
    SE_TE(nF) = SE_info{nF}.TE; % msec
end

%% ---- Get data info of ESSE----
ESSE_TE = zeros(length(ESSE_info), 1,'double');
for nF = 1:length(ESSE_info)
    ESSE_TE(nF) = ESSE_info{nF}.TE; % msec
end
t_shift = [-12 -5 -2.3  -1 0 1 2.3 5 12].' + ESSE_TE;
zero_shift_ESSE = 48;

%% Concatnate SE & ESSE
SE_idx = [2:5];
TE_all = [SE_TE(SE_idx);ESSE_TE(1)]; % all spin-echo time
echoTime_all = [SE_TE(SE_idx);t_shift]; % all echo time
cc_img_List = {SE_ccImg{SE_idx},ESSE_ccImg{1:end}};

%% Show coil images
fSIZE = 16;
XTICK = [ size(cc_img_List{1},2) * [0.5:1:size(cc_img_List,2)] ];

[ha, pos] = tight_subplot(4,1,[.01 .01],[.04 .01],[.02 .03]);
axes(ha(1));  imagesc(1e2* abs(cell2mat(cc_img_List)), [0 1]); colormap(gca,gray); axis image; colorbar; ylabel('Magnitude')
axes(ha(2));  imagesc(angle(cell2mat(cc_img_List))/pi*180,[-180 180]); colormap(gca,hsv); axis image; cb = colorbar; ylabel('Phase');ylabel(cb,'degree','FontWeight','BOLD')
axes(ha(3));  imagesc(1e2* real(cell2mat(cc_img_List)), [-1 1]); colormap(gca,gray); axis image; colorbar; ylabel('Real')
axes(ha(4));  imagesc(1e2* imag(cell2mat(cc_img_List)), [-1 1]); colormap(gca,gray); axis image; colorbar; ylabel('Imaginary')

set(ha,'XTickLabel','','FontSize',fSIZE,'FontWeight','BOLD'); set(ha,'YTickLabel','')
set(ha(4),'XTick',XTICK, 'XTickLabel',echoTime_all)


%% Complex Non-linear Least-square fitting 
load('F:\USC\MREL\LowField\LungImaging\T2measurement\Data\0528AppleGroundTruth\MAT\bw_img.mat');
temp = cc_img_List;
for k = 1: length(temp)
    signal_TE(:,k) = 1e2*temp{k}(bw_img); % each columne for one TE; each row for one voxel in ROI
end

tic
parfor k = 1:size(signal_TE,1)
    k
    y = [real(signal_TE(k,:)).'; imag(signal_TE(k,:)).'];
    x0 = [real(signal_TE(k,1)), imag(signal_TE(k,1)), 1/70, 1/13, 50].'; % inital values
    E = @(x) mixed_rho_T2_T2prime_B0_fitting(x, echoTime_all, y, TE_all, zero_shift_ESSE);
    options = optimoptions('lsqnonlin');
    [x(k,:),resnorm,residual,exitflag,output] = lsqnonlin(E,x0,[],[],options);
end
toc

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
    plot(TE, reshape(angle(im_tse(row,col,:,:))/pi*180, [length(TE) 1]),'o-');
%     xline(TE(zero_shift_ind),'r--')
%     yline(abs(im_tse(row,col,:,zero_shift_ind)),'r--')
    ylabel('Phase(degree)','FontWeight','BOLD')
    xlabel('Echo-shifted Time (ms)','FontWeight','BOLD')
    title('Phase of ES-SE data','FontWeight','BOLD')
    set(gca,'XTick',TE, 'XTickLabel',TE)
    if button == 3
        return;
    end
end

%% T2 map 
T2map_ESSE = zeros(size(bw_img));
T2map_ESSE(bw_img) = 1./x(:,3);

figure,
imagesc(T2map_ESSE,[0 100])
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
imagesc(T2primemap_ESSE,[0 30])
h = colorbar;
ylabel(h, 'T2'' (ms)','FontWeight','bold')
colormap(hot)
axis image off
title('T2'' map from the Mixed')

% Mat_file = fullfile(Mat_folder,'T2primemap_ESSE.mat'); % MAT name
% if ~isfile(Mat_file)
%     save (Mat_file,'T2primemap_ESSE');
% end

%% B0 map 
B0map_ESSE = zeros(size(bw_img));
B0map_ESSE(bw_img) = x(:,5);

figure,
imagesc(B0map_ESSE,[-60 60])
h = colorbar;
ylabel(h, '\Deltaf (Hz)','FontWeight','bold')
colormap(hsv)
axis image off
title('\Deltaf map from the Mixed')

% Mat_file = fullfile(Mat_folder,'B0map_ESSE.mat'); % MAT name
% if ~isfile(Mat_file)
%     save (Mat_file,'B0map_ESSE');
% end