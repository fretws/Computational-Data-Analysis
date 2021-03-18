%% SVD on tracked movement

%% Case 1

X1 = double([
    transpose(track1_1);
    transpose(track2_1);
    transpose(track3_1)
]);

%%
[U,~,~] = svd(X1, 'econ');
displacement = transpose(U(:, 1))*X1;

%% Plotting the tracked points projected onto the 1-dimensional basis formed by the main component

figure(1)
plot(displacement);
title("Extracted Linear Displacement, with Minimal Noise");
xlabel("Time in Frames");
ylabel("transpose(U(:,1))*X1");
axis([0 length(displacement) -inf inf]);
axis 'auto y'

%% Case 2

X2 = double([
    transpose(track1_2);
    transpose(track2_2);
    transpose(track3_2)
]);

%%
[U,~,~] = svd(X2, 'econ');
displacement = transpose(U(:, 1))*X2;

%% Plotting the tracked points projected onto the 1-dimensional basis formed by the main component

figure(2)
plot(displacement);
title("Extracted Linear Displacement, with Exaggerated Camera Shake");
xlabel("Time in Frames");
ylabel("transpose(U(:,1))*X2");
axis([0 length(displacement) -inf inf]);
axis 'auto y'

%% Case 3

X3 = double([
    transpose(track1_3);
    transpose(track2_3);
    transpose(track3_3)
]);

%%
[U,~,~] = svd(X3, 'econ');
displacement = transpose(U(:, 1))*X3;

%% Plotting the tracked points projected onto the 1-dimensional basis formed by the main component

figure(3)
plot(displacement);
title("Extracted Linear Displacement, with Off-Axis Displacement");
xlabel("Time in Frames");
ylabel("transpose(U(:,1))*X3");
axis([0 length(displacement) -inf inf]);
axis 'auto y'

%% Case 4

X4 = double([
    transpose(track1_4);
    transpose(track2_4);
    transpose(track3_4)
]);

%%
[U,~,~] = svd(X4, 'econ');
displacement = transpose(U(:, 1))*X4;
% Remove data from four frames that did not show a sufficient amount of
% yellow paint on the bucket for tracking the bucket's position.
displacement = rmoutliers(displacement);

%% Plotting the tracked points projected onto the 1-dimensional basis formed by the main component

figure(4)
plot(displacement);
title("Extracted Linear Displacement, with Rotation of Mass");
xlabel("Time in Frames");
ylabel("transpose(U(:,1))*X4");
axis([0 length(displacement) -inf inf]);
axis 'auto y'
