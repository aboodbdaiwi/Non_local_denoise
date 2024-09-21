
function denoise_image_all=denoise_nonlocal_mppca(KSP2,nordic,kernel)


mask_all = zeros(size(KSP2));

denoise_image_all = zeros(size(KSP2));


for slice=floor(kernel/2)+1:3:size(KSP2,3)-floor(kernel/2)

im_r= squeeze(KSP2(:,:,slice-floor(kernel/2):slice+floor(kernel/2),:));
%
%    im_r(:,:,1:34)= real(squeeze(dwi0_noisy(:,:,slice,1:34)));
% im_r(:,:,35:68) =imag(squeeze(dwi0_noisy(:,:,slice,1:34)));
%im_r2=zeros(250,270,32); im_r2(:,:,:)=im_r(11:260,6:275,:);
%im_r2=im_r;
%% denoise with proposed framework
% VST
rimavst= im_r;
x_image=rimavst;
% denoise using optimal shrinkage
for iter=1:1
% tic
y_image=x_image+0.2*(rimavst-x_image);
noisy_image=y_image;
step=1;
dir=size(im_r,4);
numm=0;
for i=1:step:size(noisy_image,1)-kernel+1
for j=1:step:size(noisy_image,2)-kernel+1
numm=numm+1;
end
end
noisy_patch2=zeros(kernel,kernel,kernel,dir,numm);
noisy_patch3=zeros(kernel,kernel,kernel,dir,numm);
numm=0;
for i=1:step:size(noisy_image,1)-kernel+1
for j=1:step:size(noisy_image,2)-kernel+1
numm=numm+1;
noisy_patch2(:,:,:,:,numm)=squeeze(noisy_image(i:i+kernel-1,j:j+kernel-1,:,:));
end
end
noisy_patch2=reshape(noisy_patch2,[kernel*kernel*kernel*dir,numm]);
numm=0;
for i=1:step:size(noisy_image,1)-kernel+1
for j=1:step:size(noisy_image,2)-kernel+1
numm=numm+1;
noisy_patch3(:,:,:,:,numm)=squeeze(nordic(i:i+kernel-1,j:j+kernel-1,slice-1:slice+1,:));
end
end
noisy_patch3=reshape(noisy_patch3,[kernel*kernel*kernel*dir,numm]);
D=pdist2(noisy_patch3',noisy_patch3');
k=140;
denoise_patch=[];
denoise_image=zeros(size(noisy_image));
mask=zeros(size(noisy_image));
numm2=0;
for i=1:step:size(noisy_image,1)-kernel+1
for j=1:step:size(noisy_image,2)-kernel+1
numm2=numm2+1;
[dis,number] = sort(D(:,numm2));
sorted=number(1:1+k);
sorted_patches=squeeze(noisy_patch2(:,sorted));
sorted_patches=reshape(sorted_patches,[kernel,kernel,kernel,dir,k+1]);
sorted_patches=permute(sorted_patches,[1 2 4 3 5]);
sorted_patches=reshape(sorted_patches,[kernel*kernel*dir,kernel*(k+1)]);
centering=0;
MM=size(sorted_patches,2);
NNN=size(sorted_patches,1);
R = min(MM, NNN);
scaling = (max(MM, NNN) - (0:R-centering-1)) / NNN;
scaling = scaling(:);
sorted_patches=sorted_patches';
[u,sigma,v]=svd(sorted_patches,'econ');
S=diag(sigma);
vals=S;
vals = (vals).^2 / NNN;
% First estimation of Sigma^2;  Eq 1 from ISMRM presentation
csum = cumsum(vals(R-centering:-1:1)); cmean = csum(R-centering:-1:1)./(R-centering:-1:1)'; sigmasq_1 = cmean./scaling;
% Second estimation of Sigma^2; Eq 2 from ISMRM presentation
gamma = (MM - (0:R-centering-1)) / NNN;
rangeMP = 4*sqrt(gamma(:));
rangeData = vals(1:R-centering) - vals(R-centering);
sigmasq_2 = rangeData./rangeMP;
th = find(sigmasq_2 < sigmasq_1, 1)
if isempty(th)
th=MM
end
keep=1:th-1;
%
% for ss=1:min(size(sigma,1),size(sigma,2))
%      w= 2.2* sqrt(140)/(sigma(ss,ss)+eps);% 0.0004 for hcp
%
%      sigma(ss,ss)=max(0, sigma(ss,ss)-w);
% end
%         temp= u(:,:) * sigma(:,:) * v(:,:)';
%         for th=1:size(sigma,1)
%             if sum(sum(sigma(1:th,1:th)))/sum(sum(sigma))>0.95
%                 keep=1:th break
%             end
%         end
% for ss=1:min(size(sigma,1),size(sigma,2))
%      w= 120* sqrt(140)/(sigma(ss,ss)+eps); sigma(ss,ss)=max(0,sigma(ss,ss)-w);
% end
%          temp= u(:,:) * sigma(:,:) * v(:,:)';
%
temp=u(:,keep) * sigma(keep,keep) * v(:,keep)';
temp=temp';
denoise_patch(:,1:kernel)=temp(:,1:kernel);
tempp=reshape(denoise_patch(:,1:kernel),[kernel,kernel,dir,kernel]);
tempp=permute(tempp,[1 2 4 3]);
denoise_image(i:i+kernel-1,j:j+kernel-1,:,:)=  denoise_image(i:i+kernel-1,j:j+kernel-1,:,:)+tempp(1:kernel,1:kernel,:,:)*1;
mask(i:i+kernel-1,j:j+kernel-1,:,:)=mask(i:i+kernel-1,j:j+kernel-1,:,:)+1;
end
end
denoise_image=denoise_image./(eps+mask);

%clear D
end
denoise_image_all(:,:,slice-floor(kernel/2):slice+floor(kernel/2),:) = denoise_image;

mask_all(:,:,slice-floor(kernel/2):slice+floor(kernel/2),:) =mask_all(:,:,slice-floor(kernel/2):slice+floor(kernel/2),:) +1;

end
denoise_image_all=denoise_image_all./(eps+mask_all);

end