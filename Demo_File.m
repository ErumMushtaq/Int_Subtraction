clear all; clc; close all;

%% Multi-channel brain data from AC LORAKS Demo set
load MPRAGE_multi_channel % Load k-space data

[nvx nvy nc] = size(kData);
st = kData;
ns = 4; % Size of sniffer Coil

%% Calibration Phase
% Interference Signals 
%wI(1,:) = fft(BPSK(127));
wI(1,:) = fft(bpsk_new([1 1 0 0 0], 128));  % bpsk_new(data_stream, fc)

% Channels H(w) & G(w)
 si = size(wI,1); % number of interferes
 wH = randn(nc, si)*0.2 + i*randn(nc, si)*0.2;
 wG = randn(ns, si)*0.2 + i*randn(ns, si)*0.2;

 wS = wG * wI(:,50);
 wR = wH * wI(:,50);

% Compute Transform
 TLeft_inverse = (inv(wS*wS')* wS * wR')';
 TPseudo_inv = (pinv(wS')*wR')';
 TPatent_T = wH * inv(wG' * wG ) * wG';
 
%% Acquisition Phase
fc = 200;
 for kx=1:1:nvx

   wI(1,:) = fft(bpsk_new([1 0 1 1 0], fc)); % change frequency on every itertion
   fc = fc+20;
   for ky=1:1:nvy
       
 %% Define Signals
 %% M(w)
 wM(1:1:nc) = fft(st(kx,ky,:));
 wR = wM' + wH * wI(:,kx);
 wS = wG * wI(:,kx);
 %% Plot signals
% figure;subplot(1,3,1);plot(abs(fft(st(128,:,1)))); axis square; title('M(f)'); % 128 TR
% subplot(1,3,2);plot(abs(wI(1,1:1:256))); axis square; title('I(f)'); % 128 TR
% subplot(1,3,3);plot(abs((fft(st(128,:,1)))+(wI(1,1:1:256)))); axis square; title('R(f)');
 

 %% Estimate Mw using transform (from calibration Phase)
 e_Left_inverse = wR - TLeft_inverse*wS;
 e_TPseudo_inv = wR - TPseudo_inv*wS;
 e_TPatent_T = wR - TPatent_T*wS;
 
 %% K space
 k_Left_inverse = ifft(e_Left_inverse');
 k_TPseudo_inv = ifft(e_TPseudo_inv');
 k_TPatent_T = ifft(e_TPatent_T');
 k_TR = ifft(wR');
 
 kspace_LI(kx,ky,:) = k_Left_inverse;
 kspace_PI(kx,ky,:) = k_TPseudo_inv;
 kspace_PT(kx,ky,:) = k_TPatent_T;
 kspace_wR(kx,ky,:) = k_TR;

 
   end 
 end
 I_int = ifft2c(kspace_wR);
 I_LI = ifft2c(kspace_LI); % nx, ny, nc, nt
 I_PI = ifft2c(kspace_PI); % nx, ny, nc, nt
 I_PT = ifft2c(kspace_PT); % nx, ny, nc, nt
 I_WI = ifft2c(st);
 
  %% Error
% error = I_wI - I_PI;
 %% Coil Combination
%  Sum of squares
%  rSoS = sqrt(sum(abs(fftshift(ifft2(ifftshift(recon)))).^2,3));
%  rSoS = sqrt(sum(abs(fftshift(ifft2(ifftshift(recon)))).^2,3));
%  rSoS = sqrt(sum(abs(fftshift(ifft2(ifftshift(recon)))).^2,3));
%  rSoS = sqrt(sum(abs(fftshift(ifft2(ifftshift(recon)))).^2,3));
 
 img_int(:,:) = senseR1(( I_int(:,:,:)), coil_sens, eye(nc));
 img_LI(:,:) = senseR1(( I_LI(:,:,:)), coil_sens, eye(nc));
 img_PI(:,:) = senseR1(( I_PI(:,:,:)), coil_sens, eye(nc));
 img_PT(:,:) = senseR1(( I_PT(:,:,:)), coil_sens, eye(nc));
 img_WI(:,:) = senseR1(( I_WI(:,:,:)), coil_sens, eye(nc));
 
 figure(1);
 subplot(1,4,2);imagesc(abs(img_int)); axis square; title('Interference'); axis off;
 subplot(1,4,1);imagesc(abs(img_WI)); axis square; title('Ground Truth'); axis off;
 subplot(1,4,3);imagesc(abs(img_PI)); axis square; title('Pseudo Inverse Solution'); axis off;
 subplot(1,4,4);imagesc(abs(img_PT)); axis square; title('Patent Solution'); axis off;
 
