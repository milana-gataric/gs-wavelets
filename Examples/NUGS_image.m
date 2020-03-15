% Reconstruction of an image; data is precomputed for 'earth.jpg'
% Uses precomputed data (available from http://www.damtp.cam.ac.uk/research/afha/code/).

clear all;
close all;

%% The data is precomputed for the following parametars

% R = 8; % maximum scale of wavelet coefficients
% a = 2; % number of vanishing moments

% D = sqrt(2)/4; % the bound for the density delta < D
% ss = 'polar'; % type of the sampling scheme

% N = 2^R; % number of reconstruction vectors
% K = N; % band-width of the samples, i.e. samples are taken from the ball of radius K

%% Change if you want image in color

color = 'no'; 

%% Load the data

load('polar_ss_256');

load('ft_sca_256DB2_polar256')

%% Load the image

im = imread('earth.jpg');
if strcmp(color,'no')
    im=double(rgb2gray(im));
end
S = length(im);

%% Load the samples; if you use different image then compute the samples

if strcmp(color,'no')
    %sts = nufft_init(2*pi*1/S*sp,[S S],[7 7],[2*S 2*S],[0 0]);
    %omega_plain = nufft((1/S)^2*im,sts);
    load('samples_earth_gray_polar256')
    omega = sqrt(mu).*omega_plain;
else
    %fc_red=double(im(:,:,1)); fc_green=double(im(:,:,2)); fc_blue=double(im(:,:,3));
    %sts = nufft_init(2*pi*1/S*sp,[S S],[7 7],[2*S 2*S],[0 0]);   
    %omega_plain_red = nufft((1/S)^2*fc_red,sts);
    %omega_plain_green = nufft((1/S)^2*fc_green,sts);
    %omega_plain_blue = nufft((1/S)^2*fc_blue,sts);
    load('samples_earth_color_polar256')
    omega_red = sqrt(mu).*omega_plain_red;
    omega_green = sqrt(mu).*omega_plain_green;
    omega_blue = sqrt(mu).*omega_plain_blue;
end

%% Initialize NUFFT

disp('Data loaded. Initializing NUFFT...')
stx = nufft_init(2*pi*1/2^R*sp(:,1),2^R,7,2*2^R,0,'minmax:tuned');
sty = nufft_init(2*pi*1/2^R*sp(:,2),2^R,7,2*2^R,0,'minmax:tuned');
st = nufft_init(2*pi*1/2^R*sp,[2^R 2^R],[7 7],[2*2^R 2*2^R],[0 0],'minmax:tuned');

%% Compute the GS function 

disp('Computing GS reconstruction...')
GS_handle = @(x,mode) CDJV_Op2_VecHandle_NU(x, mode, R, a, ft_sca_x, ft_sca_L_x, ft_sca_R_x, ft_sca_y, ft_sca_L_y, ft_sca_R_y, sp, mu, st, stx, sty);

if strcmp(color,'no')
    z = lsqr(GS_handle,omega(:),10^(-12),50);
    wcoeff = reshape(z, sqrt(length(z)), sqrt(length(z)));
   
elseif strcmp(color,'yes')
    z_red = lsqr(GS_handle,omega_red(:),10^(-12),50);
    wcoeff_red = reshape(z_red, sqrt(length(z_red)), sqrt(length(z_red)));

    z_green = lsqr(GS_handle,omega_green(:),10^(-12),50);
    wcoeff_green = reshape(z_green, sqrt(length(z_green)), sqrt(length(z_green)));

    z_blue = lsqr(GS_handle,omega_blue(:),10^(-12),50);
    wcoeff_blue = reshape(z_blue, sqrt(length(z_blue)), sqrt(length(z_blue)));
end
%% Reconstruction from wavelets

disp('Computing image from reconstructed wavelet coefficients...')
p=12;
if strcmp(color,'no')
    reconim=get2DReconstruction(wcoeff, a, p, R);
    reconim_nw=get2DReconstruction(wcoeff_nw, a, p, R);
elseif strcmp(color,'yes')
    color_rec_red=get2DReconstruction(wcoeff_red, a, p, R);
    color_rec_green=get2DReconstruction(wcoeff_green, a, p, R);
    color_rec_blue=get2DReconstruction(wcoeff_blue, a, p, R);

    reconim = zeros(2^p,2^p,3);
    reconim(:,:,1)=color_rec_red;
    reconim(:,:,2)=color_rec_green;
    reconim(:,:,3)=color_rec_blue;
    reconim = uint8(real(reconim));

end

%% Resize the original image

im_res=imresize(im,2^p/S,'bicubic');


%% Make figures

figure('Position', [100, 100, 1400, 600],'Name','Image reconstruction from radial measurements');
subplot(1,3,1); imshow(real(im),[0,255]); title('Original')
subplot(1,3,2); imshow(real(reconim),[0,255]); title('GS reconstruction with DB2')
subplot(1,3,3); imshow(imcomplement(abs(reconim-im_res)),[-255,0]); title('| Original - Reconstruction |')

if strcmp(color,'no')
    error_gs_rel=sqrt(sum(sum(abs(im_res - reconim).^2))/2^(2*p))/sqrt(sum(sum(abs(im_res).^2))/2^(2*p));
else
    error_gs_rel_r=sqrt(sum(sum(abs(im_res(:,:,1) - reconim(:,:,1)).^2))/2^(2*p))/sqrt(sum(sum(abs(im_res(:,:,1)).^2))/2^(2*p));
    error_gs_rel_g=sqrt(sum(sum(abs(im_res(:,:,2) - reconim(:,:,2)).^2))/2^(2*p))/sqrt(sum(sum(abs(im_res(:,:,2)).^2))/2^(2*p));
    error_gs_rel_b=sqrt(sum(sum(abs(im_res(:,:,3) - reconim(:,:,3)).^2))/2^(2*p))/sqrt(sum(sum(abs(im_res(:,:,3)).^2))/2^(2*p));
    error_gs_rel=(error_gs_rel_r+error_gs_rel_g+error_gs_rel_b)/3;
end

fprintf('Relative GS reconstruction error is %d \n',error_gs_rel); 



% Copyright (c) 2014. Clarice Poon and Milana Gataric
