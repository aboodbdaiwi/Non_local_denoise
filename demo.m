
load('MB2_AP_Gre_img_GRAPPA_PF_ACCsb_GhostMethod3_Avg1')
addpath ./utils
addpath ./NORDIC_Raw-main
addpath ./NIfTI_20140122
%% image_data is a complex data 
%% use NORDIC preprocessing   
image_data = IMG_GRAPPA(:,:,40:80,:);
path = 'D:\research\non_local_denoise_code\';
A = make_nii(abs(single(image_data(:,:,:,:))),4);
save_nii(A,[path,'mag.nii'])
A = make_nii(angle(single(image_data(:,:,:,:))),4);
save_nii(A,[path,'phase.nii'])
ARG.DIROUT = path;
ARG.temporal_phase = 3;
ARG.phase_filter_width = 3;
ARG.use_generic_NII_read = 1;

kernel = 3;

NORDIC_denoise_svs([path,'mag.nii'],[path,'pha.nii'],['denoised.nii'],ARG,kernel)

NORDIC_denoise_mppca([path,'mag.nii'],[path,'pha.nii'],['denoised.nii'],ARG,kernel)


%% or use another way pre_denoising_im refers to an initial denoising result from NORDIC or other methods
kernel = 3;

NIFTI_NORDIC([path,'mag.nii'],[path,'phase.nii'],['nordic.nii'],ARG,path)
pre_denoising_im = load_nii(['phase.nii']);
pre_denoising_im = pre_denoising_im.img;

denoise_image_all = denoise_nonlocal_svs(image_data,pre_denoising_im,kernel);  

denoise_image_all = denoise_nonlocal_mppca(image_data,pre_denoising_im,kernel);
