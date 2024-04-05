function [mpsnr,mssim,mfsim,ergas,sam]=msqia(imagery1, imagery2)
% Evaluates the quality assessment indices for two HSIs.
%
% Syntax:
%   [mpsnr, mssim,ergas ] = MSIQA(imagery1, imagery2)
% Input:
%   imagery1 - the reference HSI data array
%   imagery2 - the target HSI data array
%   NOTE: MSI data array  is a M*N*K array for imagery with M*N spatial
%	pixels, K bands and DYNAMIC RANGE [0,1];
[M,N,p]  = size(imagery1);
psnrvector=zeros(1,p);
for i=1:1:p
    J=255*imagery1(:,:,i);
    I=255*imagery2(:,:,i);
    psnrvector(i)=PSNR_me(J,I);
end 
mpsnr = mean(psnrvector);

SSIMvector=zeros(1,p);
for i=1:1:p
    J=255*imagery1(:,:,i);
    I=255*imagery2(:,:,i); 
    [ SSIMvector(i),~] = ssim(J,I);
end
mssim=mean(SSIMvector);

ergas = ErrRelGlobAdimSyn(255*imagery1, 255*imagery2);

FSIMvector=zeros(1,p);
for i=1:1:p
    J=255*imagery1(:,:,i);
    I=255*imagery2(:,:,i); 
    FSIMvector(i) = FeatureSIM(J,I);
end
mfsim = mean(FSIMvector);

sam = SpectAngMapper(255*imagery1, 255*imagery2);