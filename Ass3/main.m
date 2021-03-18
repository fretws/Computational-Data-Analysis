%{
02-20-2021
Shane Fretwell
AMATH 482 Assignment 3, Spring-Mass Systems
%}
%%
clear; close all;

%%
GRAPHICS = false;
w = 640;
h = 480;

%%
load('C:/Users/sear/Documents/AMATH482_data/Ass3/cam1_1.mat')
load('C:/Users/sear/Documents/AMATH482_data/Ass3/cam2_1.mat')
load('C:/Users/sear/Documents/AMATH482_data/Ass3/cam3_1.mat')

%%
if GRAPHICS
    implay(vidFrames1_1)
end
%%
if GRAPHICS
    implay(vidFrames2_1)
end
%%
if GRAPHICS
    implay(vidFrames3_1)
end

%% Show color masked video
if GRAPHICS
    figure(2)
    for k=1:T3
    [BW,maskedRGBImage] = cam3_1Mask(vidFrames3_1(:,:,:,k));
    imshow(maskedRGBImage);
    end
end

%% Tracking cam1_1 by color
T = 212;
track1_1 = trackByColor(vidFrames1_1,1,T,@cam1_1Mask);

%% Tracking cam2_1 by color
start = 10;
track2_1 = trackByColor(vidFrames2_1,start,T,@cam2_1Mask);

%% Tracking cam3_1 by color
start = 20;
track3_1 = trackByColor(vidFrames3_1,start,T,@cam3_1Mask);

%% Plotting tracked points
if GRAPHICS
    figure(1)
    for k=1:T
        scatter(track1_1(k, 1), track1_1(k, 2), 25);
        axis([0 w 0 h]);
        hold on;
        scatter(track2_1(k, 1), track2_1(k, 2), 25);
        scatter(track3_1(k, 1), track3_1(k, 2), 25);
        hold off;

        drawnow
    end
end

%% Move on to Case 2
clear vidFrames1_1 vidFrames2_1 vidFrames3_1
load('C:/Users/sear/Documents/AMATH482_data/Ass3/cam1_2.mat')
load('C:/Users/sear/Documents/AMATH482_data/Ass3/cam2_2.mat')
load('C:/Users/sear/Documents/AMATH482_data/Ass3/cam3_2.mat')

%%
if GRAPHICS
    implay(vidFrames1_2)
end
%%
if GRAPHICS
    implay(vidFrames2_2)
end
%%
if GRAPHICS
    implay(vidFrames3_2)
end

%% Show color masked video
if GRAPHICS
    figure(5)
    T = length(vidFrames3_2(1,1,1,:));
    for k=1:T
        [BW,maskedRGBImage] = cam3_2Mask(vidFrames3_2(:,:,:,k));
        imshow(maskedRGBImage);
    end
end
%% Tracking cam1_2 by color
start = 10;
T = length(vidFrames1_2(1,1,1,:)) - start;
track1_2 = trackByColor(vidFrames1_2,1,T,@cam1_1Mask);

%% Tracking cam2_2 by color
start = 1;
track2_2 = trackByColor(vidFrames2_2,start,T,@cam2_1Mask);

%% Tracking cam3_2 by color
start = 15;
track3_2 = trackByColor(vidFrames3_2,start,T,@cam3_2Mask);
%% Plotting tracked points
if GRAPHICS
    figure(1)
    for k=1:T
        scatter(track1_2(k, 1), track1_2(k, 2), 25);
        axis([0 w 0 h]);
        hold on;
        scatter(track2_2(k, 1), track2_2(k, 2), 25);
        scatter(track3_2(k, 1), track3_2(k, 2), 25);
        hold off;
        title(string(k))
        drawnow
        pause(0.2)
    end
end
%% Move on to Case 3
clear vidFrames1_2 vidFrames2_2 vidFrames3_2
load('C:/Users/sear/Documents/AMATH482_data/Ass3/cam1_3.mat')
load('C:/Users/sear/Documents/AMATH482_data/Ass3/cam2_3.mat')
load('C:/Users/sear/Documents/AMATH482_data/Ass3/cam3_3.mat')

%%
if GRAPHICS
    implay(vidFrames1_3)
end
%%
if GRAPHICS
    implay(vidFrames2_3)
end
%%
if GRAPHICS
    implay(vidFrames3_3)
end

%% Show color masked video
if GRAPHICS
    figure(5)
    T = length(vidFrames3_3(1,1,1,:));
    for k=1:T
        [BW,maskedRGBImage] = cam3_2Mask(vidFrames3_3(:,:,:,k));
        imshow(maskedRGBImage);
    end
end
%% Tracking cam1_3 by color
start = 3;
T3 = length(vidFrames3_3(1,1,1,:));
track1_3 = trackByColor(vidFrames1_3,start,T3,@cam1_1Mask);

%% Tracking cam2_3 by color
start = 32;
track2_3 = trackByColor(vidFrames2_3,start,T3,@cam2_1Mask);

%% Tracking cam3_3 by color
start = 1;
track3_3 = trackByColor(vidFrames3_3,start,T3,@cam3_2Mask);

%% Plotting tracked points
if GRAPHICS
    figure(1)
    for k=1:T
        scatter(track1_3(k, 1), track1_3(k, 2), 25);
        axis([0 w 0 h]);
        hold on;
        scatter(track2_3(k, 1), track2_3(k, 2), 25);
        scatter(track3_3(k, 1), track3_3(k, 2), 25);
        hold off;
        title(string(k))
        drawnow
        pause(0.2)
    end
end

%% Move on to Case 4
clear vidFrames1_3 vidFrames2_3 vidFrames3_3
load('C:/Users/sear/Documents/AMATH482_data/Ass3/cam1_4.mat')
T1 = length(vidFrames1_4(1,1,1,:));

load('C:/Users/sear/Documents/AMATH482_data/Ass3/cam2_4.mat')
T2 = length(vidFrames2_4(1,1,1,:));

load('C:/Users/sear/Documents/AMATH482_data/Ass3/cam3_4.mat')
T3 = length(vidFrames3_4(1,1,1,:));

%%
if GRAPHICS
    implay(vidFrames1_4)
end
%%
if GRAPHICS
    implay(vidFrames2_4)
end
%%
if GRAPHICS
    implay(vidFrames3_4)
end

%% Show color masked video
if GRAPHICS
    figure(5)
    T = length(vidFrames2_4(1,1,1,:));
    for k=1:T
        [BW,maskedRGBImage] = cam2_1Mask(vidFrames2_4(:,:,:,k));
        imshow(maskedRGBImage);
    end
end
%% Tracking cam1_4 by color
start = 7;
T = T1 - start;
track1_4 = trackByColor(vidFrames1_4,start,T,@cam1_1Mask);

%% Tracking cam2_4 by color
start = 15;
track2_4 = trackByColor(vidFrames2_4,start,T,@cam2_1Mask);

%% Tracking cam3_4 by color

start = 7;
track3_4 = trackByColor(vidFrames3_4,start,T,@cam3_2Mask);

%% Plotting tracked points
if GRAPHICS
    figure(1)
    for k=1:T
        scatter(track1_4(k, 1), track1_4(k, 2), 25);
        axis([0 w 0 h]);
        hold on;
        scatter(track2_4(k, 1), track2_4(k, 2), 25);
        scatter(track3_4(k, 1), track3_4(k, 2), 25);
        hold off;
        title(string(k))
        drawnow
        pause(0.1)
    end
end

%%
clear vidFrames1_4 vidFrames2_4 vidFrames3_4
