

%% generate clean ref data from large vomumn data 
%% use fsl to fit a tensor 



dti=load_nii('dti_tensor.nii.gz');
dti=dti.img;


bvals=textread('bvals.txt');


bvecs=textread('bvecs.txt');
bvecs=bvecs(:,1:288);



%%generate modelbasec dwi

jj=(bvals<100);
b0Image=slice_img1(:,:,:,jj);

b0_im=dwi_img(:,:,:,1);
b0Image=mean(b0Image,4);



b0_im=data_dwi(:,:,6,1);
b0Image=mean(b0Image,4);


b0Image=repmat(b0Image,[1 1 1 144]);
slice_img1(:,:,:,jj)=b0Image;


D_tensor=zeros(145,174,145,3,3);
           D_tensor(:,:,:,1,1)=dti(:,:,:,1);
           D_tensor(:,:,:,1,2)=dti(:,:,:,2);
           D_tensor(:,:,:,1,3)=dti(:,:,:,3);
           D_tensor(:,:,:,2,1)=dti(:,:,:,2);
           D_tensor(:,:,:,2,2)=dti(:,:,:,4);
           D_tensor(:,:,:,2,3)=dti(:,:,:,5);
           D_tensor(:,:,:,3,1)=dti(:,:,:,3);
           D_tensor(:,:,:,3,2)=dti(:,:,:,5);
           D_tensor(:,:,:,3,3)=dti(:,:,:,6);
    b0_im=    b0Image;   
     D_tensor=     dti;
     bvecs=bvecsnew;
     bvals=bvalsnew;
for i=1:16
    
    
    for j=1:size(b0_im,1)
        for m=1:size(b0_im,2)
            for sl=1:size(b0_im,3)
          D=squeeze(D_tensor(j,m,sl,:,:));
                eff_D=bvals(:,i)*bvecs(:,i)'*D*bvecs(:,i);
                dwi_generated2(j,m,sl,i)=b0_im(j,m,sl)*exp(-eff_D);
    
            end
        end
    end
end




%% generate mean dwi

b0_im=dwi_img(:,:,:,1);
for i=1:16
tic
    B_matrix=bvals(:,i);
   B_matrix=repmat(B_matrix,[145 174 145]);
   
     
                eff_D=B_matrix.*md;
    dwi_generated4(:,:,:,i)=b0_im.*exp(-eff_D);
    toc
       
end






%% create complex noise 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  create spatially varying noise

ref=dwi_generated2new/max(max(max(max(dwi_generated2new))));
dwi0=ref;


sim_img_size=size(dwi0,1);
sim_img_size2=size(dwi0,2);
sim_img_size3=size(dwi0,3);
sim_temp_size=size(dwi0,4);

gfactor_1D= hamming(sim_img_size);
gfactor_1D=gfactor_1D-min(gfactor_1D(:))+1;  % set it to be at least 1.


gfactor_1D2= hamming(sim_img_size2);
gfactor_1D2=gfactor_1D2-min(gfactor_1D2(:))+1;  % set it to be at least 1.

gfactor2D= gfactor_1D*gfactor_1D2';


gfactor_1D3= hamming(sim_img_size3);
gfactor_1D3=gfactor_1D3-min(gfactor_1D3(:))+1;  % set it to be at least 1.
for z=sim_img_size3:-1:1
    gfactor(:,:,z)=gfactor2D*gfactor_1D3(z);  % surely an easier way to expand the matrix
end

for i=2:8
    
    for ttt=1:6
        
noise_level =i*0.0008;%2-20%


IID_NOISE=noise_level*complex(randn([sim_img_size sim_img_size2 sim_img_size3 sim_temp_size ]),randn([sim_img_size sim_img_size2 sim_img_size3 sim_temp_size ])   );
%%create phase
for t= sim_temp_size:-1:1
spatial_noise(:,:,:,t)= IID_NOISE(:,:,:,t).*gfactor;
end

for slice=1:sim_img_size3
[phi pt] = genphi(sim_img_size,sim_img_size2, sim_temp_size);
for pp=1:1
phaseall(:,:,slice,pp)=phi(:,:,pp)*0.4;   %0.5 1 2
end
for pp=2:25
phaseall(:,:,slice,pp)=phi(:,:,pp);   %0.5 1 2
end
for pp=26:49
phaseall(:,:,slice,pp)=phi(:,:,pp)*2;   %0.5 1 2
end

end
dwi0_noisy=dwi0.*exp(-1i*phaseall);


dwi0_noisy=dwi0_noisy+spatial_noise;

save(['/home/yexinyu/umn/mghnew/',num2str(i),'/lv',num2str(i),'rep',num2str(ttt),'.mat'],'dwi0_noisy')
    end
end
