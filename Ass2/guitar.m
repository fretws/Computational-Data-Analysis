%%
clear; close all;
%%

[y, Fs] = audioread('Floyd.m4a'); % y := intensity, Fs := number of measurements per second
y = transpose(y);
% Making the length of y even makes analysis easier without changing results in any significant way
y = y(1:end - 1);
tr_floyd = length(y)/Fs; % record time in seconds

T = 1/Fs; % Seconds between samples       
L = length(y); % Number of samples

Tau = 42; % Number of time points per second
% The Gabor Transform will be centered at each of the following points
taustep = floor(linspace(1, length(y), floor(Tau * tr_floyd)));
k = Fs*(0:(L/2))/L; % Frequencies in Hz

%%
an = 61; bn = 63; % Note range as note numbers
af = floor(16.35*L/Fs*(2^(1/12))^an); % Indices of note range in frequency domain
bf = floor(16.35*L/Fs*(2^(1/12))^bn);
specf = zeros(bf - af + 1, length(taustep));
sigma = 2.5*10^3;

for tau=1:length(taustep)
    % Define the gabor filter with L-2 norm equal to one
    filter = sqrt(2/(sigma^2*pi))*exp(-0.5*(((1:L) - taustep(tau))/sigma).^2);
    tframe = y.*filter;
    
    freqsGT = fft(tframe);
    allfreq = abs(freqsGT/L);
    posfreq = allfreq(1:(L/2 + 1));
    posfreq(2:(end - 1)) = 2*posfreq(2:(end - 1));

    specf(:, tau) = posfreq(af:bf); % D2 to C#3
end
%% Plot spectrogram
figure(3)
pcolor(taustep / Fs, log(k(af:bf) / 16.35) / log(2^(1/12)), specf)
shading interp
set(gca,'ylim',[an bn],'Fontsize',16)
colormap(hot)
colorbar
xlabel('time (seconds)'), ylabel('frequency (k)')
