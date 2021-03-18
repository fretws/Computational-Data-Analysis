%{
01-07-2021
Shane Fretwell
AMATH 482 Assignment 1, Submarine Tracking
%}

%%
clear; close all; clc

load subdata.mat % Imports the data as the 262144x49 (space by time) matrix called subdata

L = 10; % spatial domain
n = 64; % Fourier modes
T = 49; % Time points
k = (2*pi/(2*L))*[0:(n/2 - 1) -n/2:-1]; ks = fftshift(k);

Un(:,:,:,:)=reshape(subdata,n,n,n,T);

%% Time-Intensity signals, averaged over realizations from each point in space
test=reshape(Un(1,1,1,:),[1,T]);
unitfreq=fft(test,n);

figure(1)
plot(ks, fftshift(abs(unitfreq))/max(abs(unitfreq)),'r','Linewidth',2);
set(gca,'Fontsize',16)
title('Frequency Spectrum of (1, 1, 1)')
xlabel('frequency (k)')
ylabel('|unitfreq|')
grid on, axis square

ave = zeros(1, n);
for i=1:n
   for j=1:n
      for k=1:n
          fspace = fft(reshape(Un(i,j,k,:), [1,T]), n);
          ave = ave + fspace;
      end
   end
end
ave = abs(fftshift(ave))./(n^3);
%%
figure(2)
hold off
plot(ks, ave/max(ave),'r','Linewidth',2);
hold on
plot(ks, exp(-0.3*(ks + 0.94).^2),'b','Linewidth',1);
plot(ks, 0.75*exp(-0.5*(ks + 7.54).^2),'b','Linewidth',1);
set(gca,'Fontsize',16)
title('Spectrum Averaged Over Space')
xlabel('frequency (k)')
ylabel('|mean freq strength|')
grid on, axis square

