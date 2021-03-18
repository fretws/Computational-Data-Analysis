%% 
%{
03-13-2021
Shane Fretwell
AMATH 482 Assignment 5, Foreground Extraction
%}
%%
clear; close all;

%% Read in video data
rdr = VideoReader('C:/Users/sear/Documents/AMATH482_data/Ass5/ski_drop_low.mp4');
N = rdr.NumFrames;
frames = read(rdr,[1 N]);
dt = rdr.Duration / rdr.NumFrames;
t = (0:N-1)*dt;

%% Convert to Grayscale
for i=1:length(t)
    frames(:,:,1,i) = rgb2gray(frames(:,:,:,i));
end
frames = reshape(frames(:,:,1,:),540,960,length(t));

%% Reshape to be one frame per column
X = double(reshape(frames, [], N));

%% The DMD matrices
% To save space, we will not store X1 and X2, but note the following:
% X1 = X(:,1:end-1); X2 = X(:,2:end);

%% SVD of X1 and Computation of ~S
[U, Sigma, V] = svd(X(:,1:end-1),'econ');
sigma = diag(Sigma);
clear Sigma

%% Plot singular value spectrum
figure(1)
plot(log10(sigma))
title('Singular Value Spectrum of X_1, Ski Drop')
xlabel('Singular Value Index')
ylabel('Log10 of Singular Value')

%% Truncate Rank
r = 25;
U = U(:,1:r); sigma = sigma(1:r); V = V(:,1:r);

%% Obtaining Eigenfunctions
S = U'*X(:,2:end)*V*diag(1./sigma);
[eV, D] = eig(S); % compute eigenvalues + eigenvectors
mu = diag(D); % extract eigenvalues
omega = log(mu)/dt;
Phi = U*eV;

%% Solve Initial Conditions and Generate DMD Modes
%y0 = Phi\X1(:,1);
y0 = Phi\X(:,1); % pseudoinverse to get initial conditions
u_modes = zeros(length(y0),length(t));
for iter = 1:length(t)
    u_modes(:,iter) = y0.*exp(omega*t(iter));
end

%% Plot omega
figure(2)
plot(omega, 'o')
axis([min(real(omega))-3 3 -inf inf])
hold on
xline(0, '--'); yline(0, '--')
hold off
title('Values of Omega, Ski Drop')
xlabel('Re(omega)'); ylabel('Im(omega)')

%% Separate Background and Foreground
[~,i] = min(omega);
back = real(Phi(:,i)*u_modes(i,:));
fore = X - back;
% Force all positive pixel values
R = min(fore, [], 'all');
fore = fore - R;
% Force all pixel values between 0 and 255
M = max(fore, [], 'all');
XM = max(X, [], 'all');
if (M > 255)
    back = back / M * XM;
    fore = fore / M * XM;
end
back = reshape(uint8(abs(back)), 540, 960, []);
fore = reshape(uint8(abs(fore)), 540, 960, []);

%% Show the full length video, either the extracted foreground, background, or original
figure(3)
for i=1:length(t)
    %imshow(fore(:,:,i))
    imshow(back(:,:,i))
    %imshow(frames(:,:,i))
    title(strcat('Frame ', num2str(i)))
end

%% Show several original frames and separated foregrounds side by side
figure(4)
delta = 15;
start = 400;
framenums = start+(1:3)*delta;
width = length(framenums)+1;

subplot(2, width, 1)
imshow(back(:,:,framenums(1)))
title(strcat('Background Frame ', num2str(framenums(1))))

subplot(2, width, width+1)
imshow(frames(:,:,framenums(1)))
title(strcat('Frame ', num2str(framenums(1))))

for i=1:length(framenums)
    subplot(2, width, i+1)
    imshow(fore(:,:,framenums(i)))
    title(strcat('Foreground Frame ', num2str(framenums(i))))
    
    subplot(2, width, width+i+1)
    imshow(frames(:,:,framenums(i)))
    title(strcat('Frame ', num2str(framenums(i))))
end