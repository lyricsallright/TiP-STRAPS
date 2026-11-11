clear;
clc;

%% Paths and prepared input data
codeFolder = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(codeFolder, 'functions')));
% kdata_Img: [Nx, Ny, Nc, Nt], Cartesian imaging data
% mask_Img:  [Nx, Ny, Nt], actual sampling mask
% kdata_Nav: [Nt, Nx*Nc], navigating data
% sensMap_win: [Nx, Ny, Nc], precomputed coil sensitivity maps of the
% reconstruction window
% Input arrays are already reordered and share the same coil basis.
load(fullfile(codeFolder, 'data', 'prepared_data.mat'), ...
    'kdata_Img', 'mask_Img', 'kdata_Nav', 'sensMap_win');
[Nx, Ny, Nc, Nfrms_all] = size(kdata_Img);

%% adjust for the acquisition
Nfrms = 64;
FS = 6;
overlap_frms = Nfrms - FS;
first_frm = 1;
init_groups = 1;
r = 3;
norm_abs_max = 1;

%% Trajectory mask example [Ny, Nx]
% 0: trajectory region; 1: surrounding region constrained by the SC operator.
% This centered rectangle is an illustrative mask, not an experimental ROI.
tip_mask = ones(Ny, Nx);
roi_rows = max(1, round(0.45*Ny)):min(Ny, round(0.55*Ny));
roi_cols = max(1, round(0.25*Nx)):min(Nx, round(0.75*Nx));
tip_mask(roi_rows, roi_cols) = 0;

%% Temporal constraint
temporal_hold_frames = overlap_frms;
reconOPT = struct();
reconOPT.TVcoef = 1;
reconOPT.maxiter = 5;
reconOPT.tol = 1e-5;
algoFunc = @func_PSRT_recon_DtF;

%% Sliding-window reconstruction
group_all = max(0, floor((Nfrms_all-first_frm+1-Nfrms)/FS)+1);
IMG = zeros(Nx,Ny,Nfrms_all);
time_iter = zeros(group_all,1);
for i = 1:group_all
    start_frm_idx = FS*(i-1)+first_frm;
    end_frm_idx = start_frm_idx+Nfrms-1;
    kdata_Img_win = kdata_Img(:,:,:,start_frm_idx:end_frm_idx);
    mask_Img_win = mask_Img(:,:,start_frm_idx:end_frm_idx);
    kdata_Nav_win = kdata_Nav(start_frm_idx:end_frm_idx,:);


    %%  reconstruction
    if i <= init_groups
        [vt_r_win, svd_sigma_win, ~] = svds(kdata_Nav_win, r);
        V_win = vt_r_win.';
        tic
        reconSt = algoFunc(...
            kdata_Img_win, mask_Img_win, sensMap_win, V_win, reconOPT);
        toc
        time_iter(i)=reconSt.time_all;
    else
        if temporal_hold_frames > 0
            kdata_Nav_win(1:temporal_hold_frames, :) = repmat(kdata_Nav_win(temporal_hold_frames, :), temporal_hold_frames, 1);
        end
        [vt_r_win, svd_sigma_win, ~] = svds(kdata_Nav_win, r);
        V_win = vt_r_win.';

        tic
        U0 = reconSt.U0;
        reconSt = algoFunc(...
            kdata_Img_win, mask_Img_win, sensMap_win, V_win, reconOPT, tip_mask, U0);
        toc
        time_iter(i)=reconSt.time_all;

    end

%%  display
    img = funcPSRT_reconSt2img(reconSt);
    img = cutHead(img, 1);
    if i==1
        IMG(:,:,start_frm_idx:start_frm_idx+Nfrms-1)=img;
    else
        IMG(:,:,(end_frm_idx-FS+1):end_frm_idx) = img(:,:,(end-FS+1):end);
    end
end


function outIm = cutHead(inIm, thresh)
    maxVal = max(inIm(:));
    outIm = inIm / (thresh * maxVal);
    outIm(outIm > 1) = 1;
end


