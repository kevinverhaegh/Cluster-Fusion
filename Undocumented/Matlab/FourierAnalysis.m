function [f,fourierspec]=FourierAnalysis(t,L,y)
T=t(2)-t(1);
fsamp=1/T;
NFFT=2^nextpow2(L);
fourier=fft(y,NFFT)/L;
f=fsamp/2*linspace(0,1,NFFT/2+1);
fourierspec = 2*abs(fourier(1:NFFT/2+1));