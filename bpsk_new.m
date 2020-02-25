function signal = bpsk_new(d, fc)

b = 2*d-1; % Convert unipolar data stream to bipolar data stream

T = 0.00004; % Bit duration
Eb = 1/2; % Vary this to change amplitude of the waveforms
t = linspace(0,5*T,300); % discrete time sequence between 0 and 5*T (1000 samples)
N = length(t); % Number of samples
Nsb = N/length(d); % Number of samples per bit
dd = repmat(d',1,Nsb); % replicate each bit Nsb times
bb = repmat(b',1,Nsb); dw=dd'; % Transpose the rows and columns
bw = bb';
bw = bw(:)'; % Data sequence samples
%dw = dw(:)'; 
w = sqrt(2*Eb/T)*cos(2*pi*fc*t); % carrier waveform
signal = bw.*w; % modulated waveform