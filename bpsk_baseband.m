function output = bpsk_baseband(N, T, A)
%% N length of sequence being generated
%% T bit duration
%% Input signal
Baseband_signal = double((randn(N,1)>0));
T = 0.000004; % Bit duration 4 micro second
Nsb = 1/T;  % Create Sequence of specific length corresponding to Ts.
for kk = 1:1:length(Baseband_signal)
if kk ==1
signal = (repmat(Baseband_signal(kk),1,Nsb)); % replicate each bit Nsb times
else
signal = [signal (repmat(Baseband_signal(kk),1,Nsb))];
end
end
%% Sampling of the rectangular pulse
sampled_signal = signal(1:1/T:end);

%% Signal Spectrum
Spectrum = (fft2c(sampled_signal));
PSD = (fft2c(sampled_signal)).^2;
%% Plots 
% Baseband Signal
% Factor = N*T/length(signal);
% figure;plot([Factor*1e+6:Factor*1e+6:(N)*T*1e+6],signal);
% xlabel('t (\mus)'); ylabel('Amplitude'); title('I(t)');
% ylim([-1.5 1.5]); xlim([0 (N)*T*1e+6]); hold on;
% 
% % Sampled Signal
% scatter([T/2*1e+6:1*T*1e+6:N*T*1e+6],sampled_signal, '*');
% xlabel('t (\mus)'); ylabel('Amplitude'); title('I(t)');
% ylim([-1.5 1.5]); xlim([0 (N)*T*1e+6]); legend('Baseband Signal','Sampled Signal');
% 
% % I(f)
% figure;plot([(-N/(2*T)+1/(4*T))*1e-6:(1*1e-6)/T:(96/(2*T)-1/(4*T))*1e-6],((((abs(Spectrum))))));
% xlabel('f (MHz)');
% ylabel('|I(f)|');title('I(f)');
% 
% % PSD
% figure;plot([(-N/(2*T)+1/(4*T))*1e-6:(1*1e-6)/T:(96/(2*T)-1/(4*T))*1e-6],((10*log10((abs(PSD))))));
% xlabel('f (MHz)'); ylabel('|PSD| (dB)'); %title('10log10|I(f).^2|');
output = A*(Spectrum);
end



%% Extras
%figure;plot([16.0000e-12:16.0000e-12:96*0.000004],signal);
%figure;plot([16.0000e-12*1e+6:16.0000e-12*1e+6:(N)*T*1e+6],signal);
% xtickformat('%.1f'); 
% ax = gca;
% ax.XAxis.Exponent = 6;
%figure;plot([-95/(2*T):1/T:96/(2*T)],(10*log10((Spectrum))));
% xlabel('f (Hz)');
% ylabel('Magnitude'); title('I(f)');
% xtickformat('%.1f'); 
% ax = gca;
% ax.XAxis.Exponent = 6;