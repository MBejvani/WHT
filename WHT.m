% Windowed Hilbert Transform (WHT)
clear,clc
load SignalsTypes
dt=Signals.SinIF2C.dt;
s0=Signals.SinIF2C.s;
t=Signals.SinIF2C.t;
f=Signals.SinIF2C.f;
cr=Signals.SinIF2C.cr;
% Henkel Mtrix Constraction
s = s0;
N = length(s0);
s = s(:)';
s = s(1:2*fix(N/2));
N = length(s);
l = fix(N/20); % selected parameter
% s = padarray(s,[0 L]);
s = [zeros(1,l) s zeros(1,l)];
L = 2*l;
w = gausswin(L,4); % selected parameter
w = w/sum(w);
xidx = (1:L)';
yidx = 0:(N-1);
h = xidx(:,ones(N,1)) + yidx(ones(L,1),:);
H = s(h); % Henkel Matrix of the signal
figure;
imagesc(H)
% Windowed Hilbert Transform (WHT)

FH = fft(sparse(diag(w)) * H); % STFT
q = zeros(1,L);
q([1 L/2+1]) = 1;
q(2:L/2) = 2;
FWHT = sparse(diag(q)) * FH;
WHT = ifft(FWHT); % WHT
imagesc(abs(FWHT))
imagesc(real(WHT))
% First-order 
% Derivative of WHT and Chirp-rate
df = 1/(L*dt);
DHT = ifft(sparse(diag(sqrt(-1)*2*pi*df*[0:fix(L/2)-1,0,-fix(L/2)+1:-1])) * FWHT); 
IF = (real(WHT).*imag(DHT)-real(DHT).*imag(WHT)) ./ (conj(WHT).*WHT + eps); % IF
CR = diff(IF);CR(L,:) = 0; % CR
% Envelope of the Chirp Rate 

FCR = fft(sparse(diag(w)) * CR);
acr = abs(ifft(sparse(diag(q)) * FCR))/max(w);
% Smoothing
w = gausswin(L,6); % selected parameter
w = w/sum(w);
IF = (w' * IF)/2/pi;
CR = (w' * CR)/dt/2/pi;
ACR = (w' * acr)/dt/2/pi;
ACR([1:l end-l:end]) = 1/dt; % selected parameter
plot(IF)
plot(CR)
plot(ACR)
% STFT Via Chirp Rate
sd = 1./(sqrt(2*pi*ACR)+eps)/dt;
S = repmat(s0',1,N);
W = WinMtx( N,sd );
TFR = abs(fft(W.*S));
figure
imagesc(W)
imagesc(t,f,fftshift(TFR,1))
xlabel Time(s)
ylabel Frequency(Hz)
% STFT Via Original Chirp Rate
acr = abs(hilbert(cr));
acr([1:24 end-24:end]) = .45*acr([1:24 end-24:end])+875;
sd = 1./(sqrt(2*pi*acr)+eps)/dt;
W = WinMtx( N,sd );
TFR0 = abs(fft(W.*S));
figure
imagesc(W)
imagesc(t,f,fftshift(TFR0,1))
xlabel Time(s)
ylabel Frequency(Hz)
Some useful functions
% [upperEnv,lowerEnv] = envelope(x, n, method);
% wp = unwrap(s,pi,1);
% aa = angle(Fas);
% pr = phase(hilbert(sr));
% hilbert(s);
%[phi,w] = phasez(biquad)