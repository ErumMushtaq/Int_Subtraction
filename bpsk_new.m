function signal = bpsk_new(d, fc)
%d=[1 0 1 1 0]; % Data sequence
b = 2*d-1; % Convert unipolar to bipolar
T = 0.00004; % Bit duration
Eb = 1/2; % This will result in unit amplitude waveforms
%fc=0.0944/T; % Carrier frequency
t = linspace(0,5,300); % discrete time sequence between 0 and 5*T (1000 samples)
N = length(t); % Number of samples
Nsb = N/length(d); % Number of samples per bit
dd = repmat(d',1,Nsb); % replicate each bit Nsb times
bb = repmat(b',1,Nsb); dw=dd'; % Transpose the rows and columns
dw = dw(:)'; 
% Convert dw to a column vector (colum by column) and convert to a row vector
bw = bb';
bw = bw(:)'; % Data sequence samples
w = sqrt(2*Eb/T)*cos(2*pi*fc*t); % carrier waveform
signal = bw.*w; % modulated waveform