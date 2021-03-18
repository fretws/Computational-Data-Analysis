clear; close all;

%%
% Clean workspace
clear; close all; clc

load ../../AMATH482_data/Ass1/subdata.mat % Imports the data as the 262144x49 (space by time) matrix called subdata

L = 10; % spatial domain
n = 64; % Fourier modes
T = 49; % Time points
k = (2*pi/(2*L))*[0:(n/2 - 1) -n/2:-1]; ks = fftshift(k);

Un(:,:,:,:)=reshape(subdata,n,n,n,T);

%% Gabor to mark possible hits (freq strength > 0.2)
% Matrix of hit markers at every position and time
hits1 = [];
hits2 = [];

hitthreshold = .46;

b = 0.2; % filter width

freqflt1 = exp(-b*(ks + 0.94).^2);
freqflt2 = exp(-b*(ks + 7.54).^2);

a = .005;
for tau=1:T
    filter = exp(-a*((1:T) - tau).^2); % Define the filter
    % Apply the filter to the signal at all points in space
    UnG=zeros(n,n,n,T);
    for t=1:T
        UnG(:,:,:,t) = Un(:,:,:,t) .* filter(t);
    end
    for i = 1:64
       for j = 1:64
          for l = 1:64
             timeseries = reshape(UnG(i, j, l, :), [1,T]);
             unitfreq = fft(timeseries,n);
             filteredfreq1 = unitfreq.*freqflt1;
             filteredfreq2 = unitfreq.*freqflt2;
    
             processed1 = abs(ifft(filteredfreq1, T));
             processed2 = abs(ifft(filteredfreq2, T));
             
             if max(processed1) > hitthreshold
                hits1(size(hits1, 1) + 1, 1:4) = [i,j,l,tau];
             end
             if max(processed2) > hitthreshold
                hits2(size(hits2, 1) + 1, 1:4) = [i,j,l,tau];
             end
          end
       end
    end
end

%% Plot hits
figure(1)
s1 = size(hits1, 1);
x1 = hits1(:, 1); y1 = hits1(:, 2); z1 = hits1(:, 3);
time1 = hits1(:, 4);
scatter3(x1, y1, z1, 10, time1, 'filled');
title('Hits Over Time at Frequency -0.94')
xlabel('x')
ylabel('y')
zlabel('z')
cb = colorbar('Ticks', 1:4:T,...
         'TickLabels', 0:2:T/2);
cb.Label.String = 'Hour';

figure(2)
s2 = size(hits2, 1);
x2 = hits2(:, 1); y2 = hits2(:, 2); z2 = hits2(:, 3);
time2 = hits2(:, 4);
scatter3(x2, y2, z2, 10, time2, 'filled');
title('Hits Over Time at Frequency -7.54')
xlabel('x')
ylabel('y')
zlabel('z')
cb = colorbar('Ticks', 1:4:T,...
         'TickLabels', 0:2:T/2);
cb.Label.String = 'Hour';