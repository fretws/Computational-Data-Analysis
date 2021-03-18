%{
02-04-2021
Shane Fretwell
AMATH 482 Assignment 2, Music Isolation
%}
%%
clear; close all;
%%

[y, Fs] = audioread('GNR.m4a'); % y := intensity, Fs := number of measurements per second
y = transpose(y);
tr_gnr = length(y)/Fs; % record time in seconds

T = 1/Fs; % Seconds between samples       
L = length(y); % Number of samples

Tau = 42; % Number of time points per second
% The Gabor Transform will be centered at each of the following points
taustep = floor(linspace(1, length(y), floor(Tau * tr_gnr)));
k = Fs*(0:(L/2))/L; % Frequencies in Hz

%% Defining note frequencies
%{
Formula for calculating these values comes from
https://pages.mtu.edu/~suits/NoteFreqCalcs.html
%}
f0 = 16.35;
% a = 2^(1/12);
% notes = f0*a.^(0:107); % Frequency in Hz of all notes from C_0 to B_8

%% Take Gabor Transform over time
% 10538 is the length required to get notes from D#3 to A5, because D#3 is
% at index 12801 in the positive frequency space, A5 at index 2263, and
% 12800-2263=10537
specf = zeros(10538, length(taustep));
sigma = 4*10^3;
a = 5*10^-5;

for tau=1:length(taustep)
    % Define the gabor filter with L-2 norm equal to one
    filter = sqrt(2/(sigma^2*pi))*exp(-0.5*(((1:L) - taustep(tau))/sigma).^2);
    tframe = y.*filter;
    
    freqsGT = fft(tframe);
    allfreq = abs(freqsGT/L);
    posfreq = allfreq(1:(L/2 + 1));
    posfreq(2:(end - 1)) = 2*posfreq(2:(end - 1));

    specf(:, tau) = posfreq(2263:12800); % D#3 to A5
end

%% Plot spectrogram
figure(3)
pcolor(taustep / Fs, log(k(2263:12800) / f0) / log(2^(1/12)), specf)
shading interp
set(gca,'ylim',[40 70],'Fontsize',16)
colormap(hot)
colorbar
xlabel('time (seconds)'), ylabel('frequency (k)')




