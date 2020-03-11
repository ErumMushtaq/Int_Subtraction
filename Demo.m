clear all; close all;

%% Load K_space and Coil_sensitivty Maps
load('kspace.mat');
load('cmaps.mat');
[nvx nvy nc] = size(kspace);
st = kspace;
si = 1;
ns = 4; % Size of sniffer Coil
% Channels without path loss
wH = randn(nc, si) + i*randn(nc, si);
wG = randn(nc, si) + i*randn(nc, si);
%% Plots
T = 0.000004; % Bit duration 4 micro second
MM = fft2c(kspace(:,46,1));  % OPERATE IN FREQUENCY DOMAIN
wI(1,:) = bpsk_baseband(nvx, T);
% % plot([T/2:1*T:96*T],abs((kspace(:,:,1))));
% % xlabel('t (sec)');
% % ylabel('Magnitude'); title('s(t)'); grid on;
% figure; 
% % %plot(abs(MM)/length(MM));
% 
% 
% 
% plot([-96/(2*T)+1/(4*T):1/T:96/(2*T)-1/(4*T)],abs(MM));
% xlabel('f (Hz)');
% ylabel('Magnitude'); 
% title('S(f)'); grid on;
% % xtickformat('%.1f') 
% % ax = gca; ax.XAxis.Exponent = 6;
% hold on;
% plot([-96/(2*T)+1/(4*T):1/T:96/(2*T)-1/(4*T)],abs(wI));
% xlabel('f (Hz)');
% ylabel('Magnitude');
% hold on;
% wR = transpose(MM) + (wI);
% plot([-96/(2*T)+1/(4*T):1/T:96/(2*T)-1/(4*T)],abs(wR));
% xlabel('f (Hz)');
% ylabel('Magnitude'); title('S(f) + I(f)'); grid on;
coil_sens = FS_cmaps;
%% Calibration Phase
 wI(1,:) = bpsk_baseband(nvx, T);
% Channels H(w) & G(w)
 si = size(wI,1); % number of interferes
%  wH = randn(nc, si)*0.1 + i*randn(nc, si)*0.1;
%  wG = randn(ns, si)*0.1 + i*randn(ns, si)*0.1;
 wS = wG * wI(:,50);
 wR = wH * wI(:,50);

% Compute Transform
% TLeft_inverse = (inv(wS*wS')* wS * wR')';
 TPseudo_inv = (pinv(wS')*wR')';
% TPatent_T = wH * inv(wG' * wG ) * wG';

 M(:,:,:) = fft2c(st(:,:,:));  % Operate in frequency domain
 for kk = 1:1:nvx
     wI(kk,:) = (bpsk_baseband(nvx, T)); % change frequency on every itertion
 end

%% Acquisition Phase
 for kx=1:1:nvx
   for ky=1:1:nvy      
 %% Define Signals
 %% M(w)
 wM(:) = M(kx, ky,:);
% wR = transpose(wM) +  (repmat(wI(kx,ky),4,1));
 wR = transpose(wM) +  wH * (wI(kx,ky));
 wS = wG * wI(kx,ky);
% % e_Left_inverse(kx, ky, :) = wR - TLeft_inverse*wS;
e_TPseudo_inv(kx,ky,:) = wR - TPseudo_inv*wS;
% e_TPatent_T(kx,ky,:) = wR - TPatent_T*wS;
e_Tnew(kx,ky, :) = wR;
   end
 end
 img_org(:,:) = senseR1((ifft2c(ifft2c(M(:,:,:)))), coil_sens, eye(nc));
 img_int(:,:) = senseR1((ifft2c(ifft2c(e_Tnew(:,:,:)))), coil_sens, eye(nc));
 img_est(:,:) = senseR1((ifft2c(ifft2c(e_TPseudo_inv(:,:,:)))), coil_sens, eye(nc));
 figure;
 subplot(3,8,1);imagesc(abs(coil_sens(:,:,1))); axis square; title('Coil Sensitivity Map 1'); axis off;colormap gray;
 subplot(3,8,2);imagesc(abs(coil_sens(:,:,2))); axis square; title('Coil Sensitivity Map 2'); axis off;colormap gray;
 subplot(3,8,3);imagesc(abs(coil_sens(:,:,3))); axis square; title('Coil Sensitivity Map 3'); axis off;colormap gray;
 subplot(3,8,4);imagesc(abs(coil_sens(:,:,4))); axis square; title('Coil Sensitivity Map 4'); axis off;colormap gray;
 subplot(3,8,5);imagesc(abs(coil_sens(:,:,5))); axis square; title('Coil Sensitivity Map 5'); axis off;colormap gray;
 subplot(3,8,6);imagesc(abs(coil_sens(:,:,6))); axis square; title('Coil Sensitivity Map 6'); axis off;colormap gray;
 subplot(3,8,7);imagesc(abs(coil_sens(:,:,7))); axis square; title('Coil Sensitivity Map 7'); axis off;colormap gray;
 subplot(3,8,8);imagesc(abs(coil_sens(:,:,8))); axis square; title('Coil Sensitivity Map 8'); axis off;colormap gray;
 
 subplot(3,8,9);imagesc(abs(ifft2c(ifft2c(e_Tnew(:,:,1))))); axis square; title('Coil 1 Image with Interference'); axis off;colormap gray;
 subplot(3,8,10);imagesc(abs(ifft2c(ifft2c(e_Tnew(:,:,2))))); axis square; title('Coil 2 Image with Interference'); axis off;colormap gray;
 subplot(3,8,11);imagesc(abs(ifft2c(ifft2c(e_Tnew(:,:,3))))); axis square; title('Coil 3 Image with Interference'); axis off;colormap gray;
 subplot(3,8,12);imagesc(abs(ifft2c(ifft2c(e_Tnew(:,:,4))))); axis square; title('Coil 4 Image with Interference'); axis off;colormap gray;
 subplot(3,8,13);imagesc(abs(ifft2c(ifft2c(e_Tnew(:,:,5))))); axis square; title('Coil 5 Image with Interference'); axis off;colormap gray;
 subplot(3,8,14);imagesc(abs(ifft2c(ifft2c(e_Tnew(:,:,6))))); axis square; title('Coil 6 Image with Interference'); axis off;colormap gray;
 subplot(3,8,15);imagesc(abs(ifft2c(ifft2c(e_Tnew(:,:,7))))); axis square; title('Coil 7 Image with Interference'); axis off;colormap gray;
 subplot(3,8,16);imagesc(abs(ifft2c(ifft2c(e_Tnew(:,:,8))))); axis square; title('Coil 8 Image with Interference'); axis off;colormap gray;
 

 subplot(3,8,17);imagesc(abs(ifft2c(ifft2c(e_TPseudo_inv(:,:,1))))); axis square; title('Interference Cancelled Coil 1'); axis off;colormap gray;
 subplot(3,8,18);imagesc(abs(ifft2c(ifft2c(e_TPseudo_inv(:,:,2))))); axis square; title('Interference Cancelled Coil 2'); axis off;colormap gray;
 subplot(3,8,19);imagesc(abs(ifft2c(ifft2c(e_TPseudo_inv(:,:,3))))); axis square; title('Interference Cancelled Coil 3'); axis off;colormap gray;
 subplot(3,8,20);imagesc(abs(ifft2c(ifft2c(e_TPseudo_inv(:,:,4))))); axis square; title('Interference Cancelled Coil 4'); axis off;colormap gray;
 subplot(3,8,21);imagesc(abs(ifft2c(ifft2c(e_TPseudo_inv(:,:,5))))); axis square; title('Interference Cancelled Coil 5'); axis off;colormap gray;
 subplot(3,8,22);imagesc(abs(ifft2c(ifft2c(e_TPseudo_inv(:,:,6))))); axis square; title('Interference Cancelled Coil 6'); axis off;colormap gray;
 subplot(3,8,23);imagesc(abs(ifft2c(ifft2c(e_TPseudo_inv(:,:,7))))); axis square; title('Interference Cancelled Coil 7'); axis off;colormap gray;
 subplot(3,8,24);imagesc(abs(ifft2c(ifft2c(e_TPseudo_inv(:,:,8))))); axis square; title('Interference Cancelled Coil 8'); axis off;colormap gray;
 
figure;
 subplot(2,3,1);imagesc(abs((img_org(:,:)))); axis square; title('Interference Off'); axis off;
 subplot(2,3,2);imagesc(abs((img_int(:,:)))); axis square; title('Interference On (Coil Combined)'); axis off;colormap gray;
 subplot(2,3,3);imagesc(abs((img_org(:,:)-img_int(:,:)))); axis square; title('Difference'); axis off;
 subplot(2,3,4);imagesc(abs((img_org(:,:)))); axis square; title('Interference Off'); axis off;
 subplot(2,3,5);imagesc(abs((img_est(:,:)))); axis square; title('Estimated Image (Coil Combined)'); axis off;colormap gray;
 subplot(2,3,6);imagesc(abs((img_org(:,:)-img_est(:,:)))); axis square; title('Difference'); axis off;
 figure;
 Im = ifft2c(ifft2c(M));
 subplot(1,5,1);imagesc(abs(Im(:,:,1))); axis square; title('Coil 1 Image'); axis off;colormap gray;
 subplot(1,5,2);imagesc(abs(Im(:,:,2))); axis square; title('Coil 2 Image'); axis off;colormap gray;
 subplot(1,5,3);imagesc(abs(Im(:,:,3))); axis square; title('Coil 3 Image'); axis off;colormap gray;
 subplot(1,5,4);imagesc(abs(Im(:,:,1))); axis square; title('Coil 4 Image'); axis off;colormap gray;
 subplot(1,5,5);imagesc(abs(img_org(:,:,1))); axis square; title('Coil 1 Image'); axis off;colormap gray;



