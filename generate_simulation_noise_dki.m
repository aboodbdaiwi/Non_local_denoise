addpath ./utils

%% generate clean ref data from large vomumn data 

%% DKI
addpath(genpath('./fanDTasia'))

load('slice84_103_all.mat')
          
bvals=textread('sessionall.bval');

bvecs=textread('sessionall.bvec');
          
mask=load_nii('k_brain_mask.nii.gz');   
mask=mask.img;
 index=(bvals==0);
  index2=(bvals==1000);
 index3=(bvals==2500);
 
 %%   method1
dwi=slice_img1(41:end-40,26:end-25,:,:);
clear slice_img1
S0=double(dwi(41:end-40,26:end-25,:,index));

S0=mean(S0,4);
 S_1real=double(dwi(41:end-40,26:end-25,:,index2));
 S_2real=double(dwi(41:end-40,26:end-25,:,index3));
 clear dwi
 GradientOrientations_1=bvecs(1:3,index2);
 BValue_1=1000;
 GradientOrientations_2=bvecs(1:3,index3);
 BValue_2=2500;
[DKI_D,DKI_W]=DEMO_DKI_Estimation_Method3_in(S0,S_1real,S_2real,GradientOrientations_1',BValue_1,GradientOrientations_2',BValue_2);


     for n=1:size(S0,1)
         tic
        for m=1:size(S0,2)
            for sl=1:20
                        dki(1:6,1)=    DKI_D(:,n,m,sl); %The 6 unique coefficients of the diffusion tensor D
     dki(7:21,1)=    DKI_W(:,n,m,sl); %The 6 unique coefficients of the diffusion tensor D
      %    dki(7:21,1)=    0; %The 6 unique coefficients of the diffusion tensor D

                dwi_generated2(n,m,sl,1)=abs(S0(n,m,sl)*exp(0));
                              dwi_generated2(n,m,sl,2:1261)=abs(S0(n,m,sl)*exp(Gbig*dki));

            end
        end
        toc
     end
     

           dwi_generated2new(:,:,:,1)=dwi_generated2(:,:,:,1);
                      dwi_generated2new(:,:,:,2:25)=dwi_generated2(:,:,:,2:25);
           dwi_generated2new(:,:,:,26:49)=dwi_generated2(:,:,:,422:445);


 bvalsnew(:,1:1)=bvals(:,1);
  bvalsnew(:,2:421)=bvals(:,index2);
 bvalsnew(:,422:1261)=bvals(:,index3);


 bvecsnew(:,1:1)=bvecs(:,1);
  bvecsnew(:,2:421)=bvecs(:,index2);
 bvecsnew(:,422:1261)=bvecs(:,index3); 


%% method2
 temp(:,:,:,1)=S0;
  temp(:,:,:,2:421)=S_1real;
 temp(:,:,:,422:1261)=S_2real;
clear S0 S_1real S_2real
 bvalsnew(:,1:1)=bvals(:,1);
  bvalsnew(:,2:421)=bvals(:,index2);
 bvalsnew(:,422:1261)=bvals(:,index3);


 bvecsnew(:,1:1)=bvecs(:,1);
  bvecsnew(:,2:421)=bvecs(:,index2);
 bvecsnew(:,422:1261)=bvecs(:,index3); 
 
 grad(:,1:3)=bvecsnew';
grad(:,4)=bvalsnew';

RobustDKIFitting(temp, grad, mask);


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
