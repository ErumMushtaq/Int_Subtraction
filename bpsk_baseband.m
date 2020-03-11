function output = bpsk_baseband(N, T )
%% N length of sequence being generated
%% T bit duration
%% Input signal
Baseband_signal = (randn(96,1)>0);
T = 0.000004; % Bit duration 4 micro second
Nsb = 1/T;  % Create Sequence of specific length corresponding T length duration.
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
Spectrum = fft2c(sampled_signal);

%% Plots 
% figure;plot([16.0000e-12:16.0000e-12:96*0.000004],signal);
% xlabel('t (sec)');
% ylabel('Magnitude'); title('I(t) Before Sampling');
% ylim([-2 2]);
% hold on;
% scatter([T/2:1*T:96*T],sampled_signal, '*');
% xlabel('t (sec)');
% ylabel('Magnitude'); title('I(t)');
% ylim([-0.5 1.5]);
% xlim([0 96*T]);
% legend('Baseband Signal','Sampled Signal');
% figure;plot([-95/(2*T):1/T:96/(2*T)],(abs(Spectrum))/length(Spectrum));
% xlabel('f (Hz)');
% ylabel('Magnitude'); title('I(f)');
% xtickformat('%.1f') 
% ax = gca;
% ax.XAxis.Exponent = 6;
output = 100*(Spectrum);
end