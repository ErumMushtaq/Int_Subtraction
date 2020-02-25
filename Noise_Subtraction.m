clear all; clc; close all;

%% Multi-channel brain data from AC LORAKS Demo set
load MPRAGE_multi_channel % Load k-space data

[nvx nvy nc] = size(kData);
st = kData;
ns = 4; % Size of sniffer Coil

%% Calibration Phase
% Interference Signals 
%wI(1,:) = fft(BPSK(127));
wI(1,:) = fft(bpsk_new([1 0 1 1 0], 128));
% wI(1,:) = 17000;
% Channels H(w) & G(w)
 si = size(wI,1); 
 wH = randn(nc, si)*0.2 + i*randn(nc, si)*0.2;
 wG = randn(ns, si)*0.2 + i*randn(ns, si)*0.2;
 
%  wS = wG * wI(:,7);
%  wR = wH * wI(:,7);

 wS = wG * 16000;
 wR = wH * 16000;

% Transform
 TLeft_inverse = (inv(wS*wS')* wS * wR')';
 TPseudo_inv = (pinv(wS')*wR')';
 TPatent_T = wH * inv(wG' * wG ) * wG';
fc = 200;
 for kx=1:1:nvx

   wI(1,:) = fft(bpsk_new([1 0 1 1 0], 700));
   fc = fc+100;
   for ky=1:1:nvy
 %% Define Signals
 %% M(w)
 wM(1:1:nc) = fft(st(kx,ky,:));
 
 %% Interference Signals 

% wI(1,:) = fft(BPSK(20));
%  plot(abs(wI(1,:)))
% wI(2,:) = fft(bpsk([1 0 1 1 0], 128));
% wI(3,:) = fft(bpsk([1 0 1 1 0], 126.5));
% wI(4,:) = fft(bpsk([1 0 1 1 0], 127.5));
% wI(5,:) = fft(bpsk([1 0 1 1 0], 125.5));
% wI(6,:) = fft(bpsk([1 0 1 1 0], 128.5));
 
 %% H(w) & G(w)
 si = size(wI,1); 
%  wH = randn(nc, si)*0.2 + i*randn(nc, si)*0.2;
%  wG = randn(ns, si)*0.2 + i*randn(ns, si)*0.2;
 
  wR = wM' + wH * wI(:,kx);
% if ky > 100 && ky << 150
%  wR = wM' + wH * wI(:,kx);
% else 
%  wR = wM';
% end
% figure;plot(abs(fft(st(:,128,1))));
% hold on;
% plot(abs(wI(1,1:256)));
% figure;subplot(1,3,1);plot(abs(wM'));
% subplot(1,3,2);plot(abs(wI));
% subplot(1,3,3);plot(abs(wR));
 wS = wG * wI(:,kx);

 %% Estimate Mw
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
% wM' - wMe
% wM' - wwMe
% wM' - wwwMe 
 
   end 
 end
 I_int = ifft2c(kspace_wR);
 I_LI = ifft2c(kspace_LI); % nx, ny, nc, nt
 I_PI = ifft2c(kspace_PI); % nx, ny, nc, nt
 I_PT = ifft2c(kspace_PT); % nx, ny, nc, nt
 I_WI = ifft2c(st);
 
 %% Coil Combination
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
 