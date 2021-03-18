%% 
%{
03-06-2021
Shane Fretwell
AMATH 482 Assignment 4, Hand-Written Digit Recognition
%}
%%
clear; close all;
%%
[Simages, Slabels] = mnist_parse(...
    'C:/Users/sear/Documents/AMATH482_data/Ass4/train-images-idx3-ubyte', ...
    'C:/Users/sear/Documents/AMATH482_data/Ass4/train-labels-idx1-ubyte');
[Timages, Tlabels] = mnist_parse(...
    'C:/Users/sear/Documents/AMATH482_data/Ass4/t10k-images-idx3-ubyte', ...
    'C:/Users/sear/Documents/AMATH482_data/Ass4/t10k-labels-idx1-ubyte');

%% SVD
A = reshape(Simages, [28^2, 60000]);
X = reshape(Timages, [28^2, 10000]);
m = mean(A,2);
A = double(A) - m;
X = double(X) - mean(X,2);
[U,S,V] = svd(A, 'econ');
s = diag(S);

%% Singular Value Spectrum
figure(1)
plot(s);
title('Spectrum of Singular Values of A')
xlabel('Singular value index'); ylabel('Singular Value')

%% 
imgnum = 3;
img = reshape(uint8(A(:,imgnum)), [28, 28]);
ranks = [784, 25, 65, 105];
figure(2)

subplot(2,2,1)
imshow(img)
title("Original Image Minus Pixel Mean")

for digit=2:4
    subplot(2,2,digit)
    % LRA = Low Rank Approximation
    LRA = U(:, 1:ranks(digit))...
        * diag(s(1:ranks(digit)))...
        * transpose(V(:, 1:ranks(digit)));
    LRA = uint8(reshape(LRA,[28,28,60000]));
    imshow(LRA(:,:,imgnum))
    title(strcat("Rank-", string(ranks(digit)), " Approximation"))
end

%% Low Rank Approximation of Several Digits
imgnums = [2, 6, 3, 12, 18];
figure(3)
for digit=1:length(imgnums)
    subplot(3, length(imgnums), digit)
    imshow(reshape(uint8(A(:,imgnums(digit))), [28, 28]))
    title("Original Image Minus Pixel Mean")
    
    subplot(3, length(imgnums), digit + length(imgnums))
    LRA = U(:, 1:25) * diag(s(1:25)) * V(:, 1:25)';
    LRA = uint8(reshape(LRA,[28,28,60000]));
    imshow(LRA(:,:,imgnums(digit)))
    title("Rank-25 Approximation")
    
    subplot(3, length(imgnums), digit + 2*length(imgnums))
    LRA = U(:, 1:45) * diag(s(1:45)) * V(:, 1:45)';
    LRA = uint8(reshape(LRA,[28,28,60000]));
    imshow(LRA(:,:,imgnums(digit)))
    title("Rank-45 Approximation")
end

%% View strokes
L = 3;
W = 5;
start = 1;
strokes = start:(start + L*W - 1);
figure(4)
for j=1:L
    for digit=1:W
        subplot(L, W, digit + W*(j-1))

        strokeMask = zeros(1,size(V,2));
        strokeMask(strokes(digit + W*(j-1))) = 1;
        strokeimg = U(:, strokes(digit + W*(j-1)))...
            * s(strokes(digit + W*(j-1)));
        strokeimg = uint8(reshape(strokeimg,[28,28]));
        imshow(strokeimg)
        title(strcat("Stroke ", string(strokes(digit + W*(j-1)))))
    end 
end

%% Project Digits onto V-Modes
figure(5)
strokes = [2, 3, 4];
n = 5000;
hold off
for digit=[1,3,4]
    SdigitMask = find(Slabels(1:n) == digit);
    scatter3(s(strokes(1))*V(SdigitMask,strokes(1)),...
        s(strokes(2))*V(SdigitMask,strokes(2)),...
        s(strokes(3))*V(SdigitMask,strokes(3)), ...
         'DisplayName', sprintf('%i',digit), 'Linewidth', 2)
    xlabel(strcat('V-Mode ', string(strokes(1))))
    ylabel(strcat('V-Mode ', string(strokes(2))))
    zlabel(strcat('V-Mode ', string(strokes(3))))
    title('Projection onto V-modes 2, 3, 4')
    hold on
end
legend

%% Two Digit Classifier
SdigitMask = find((Slabels == 1) | (Slabels == 4));
TdigitMask = find((Tlabels == 1) | (Tlabels == 4));

% Test data in the Principal Component Basis of Training Data
PCBX = U'*X; % U'*A = S*V'
PCBA = U'*A;
n = 25;
strokes = 1:n;
% V-modes specified by strokes
testfeatures = transpose(PCBX(strokes,TdigitMask));
trainfeatures = transpose(PCBA(strokes,SdigitMask));

class = classify(trainfeatures, trainfeatures, Slabels(SdigitMask), 'linear');
trainaccuracy2digit = length(find(Slabels(SdigitMask) == class))/length(class)
class = classify(testfeatures, trainfeatures, Slabels(SdigitMask), 'linear');
testaccuracy2digit = length(find(Tlabels(TdigitMask) == class))/length(class)
clear trainaccuracy2digit testaccuracy2digit;

%% Three Digit Classifier
SdigitMask = find((Slabels == 1) | Slabels == 3 | (Slabels == 4));
TdigitMask = find((Tlabels == 1) | Tlabels == 3 | (Tlabels == 4));

% Test data in the Principal Component Basis of Training Data
PCBX = U'*X; % U'*A = S*V'
PCBA = U'*A;
n = 25;
strokes = 1:n;
% V-modes specified by strokes
testfeatures = transpose(PCBX(strokes,TdigitMask));
trainfeatures = transpose(PCBA(strokes,SdigitMask));

class = classify(trainfeatures, trainfeatures, Slabels(SdigitMask), 'linear');
trainaccuracy3digit = length(find(Slabels(SdigitMask) == class))/length(class)
class = classify(testfeatures, trainfeatures, Slabels(SdigitMask), 'linear');
testaccuracy3digit = length(find(Tlabels(TdigitMask) == class))/length(class)
clear trainaccuracy3digit testaccuracy3digit;

%% Find most and least difficult digits to classify
trainaccuracies = zeros(10, 10);
testaccuracies = zeros(10, 10);
% Test data in the Principal Component Basis of Training Data
PCBX = U'*X; % U'*A = S*V'
PCBA = U'*A;
n = 25;
strokes = 1:n;

for i=0:8
    for j=i+1:9
        SdigitMask = find((Slabels == i) | (Slabels == j));
        TdigitMask = find((Tlabels == i) | (Tlabels == j));

        % V-modes specified by strokes
        testfeatures = transpose(PCBX(strokes,TdigitMask));
        trainfeatures = transpose(PCBA(strokes,SdigitMask));

        class = classify(trainfeatures, trainfeatures, Slabels(SdigitMask), 'linear');
        trainaccuracies(i+1,j+1) = length(find(Slabels(SdigitMask) == class))/length(class);
        class = classify(testfeatures, trainfeatures, Slabels(SdigitMask), 'linear');
        testaccuracies(i+1,j+1) = length(find(Tlabels(TdigitMask) == class))/length(class);
    end
end
% Subtract one from i,j to get the actual digits, due to matlab indexing
[i,j]=find(testaccuracies == min(testaccuracies(testaccuracies ~= 0),[],'all'));
min_test_accuracy = [i-1,j-1]
[i,j]=find(testaccuracies == max(testaccuracies,[],'all'));
max_test_accuracy = [i-1,j-1]
clear min_test_accuracy max_test_accuracy

%% Decision Tree Analysis
n = 25;
strokes = 1:n;

PCBX = U'*X; % U'*A = S*V'
LPCBX = PCBX(strokes,:); % Projections onto V-modes specified by strokes
PCBA = U'*A;
LPCBA = PCBA(strokes,:);

maxLPCBA = max(LPCBA, [], 'all');
features = transpose(LPCBA) / maxLPCBA;
tree=fitctree(features,Slabels,'MaxNumSplits',1000);

trainClass = predict(tree, features);
testClass = predict(tree, transpose(LPCBX) / maxLPCBA);

trainAccuracy = length(find(Slabels == trainClass))/length(trainClass)
testAccuracy = length(find(Tlabels == testClass))/length(testClass)

%% Get training and test accuracy for several different ranks and max splits
PCBX = U'*X; % U'*A = S*V'
PCBA = U'*A;
N = [4, 25, 50, 80, 100];
B = [100, 500, 1000, 1500];
treeAccuracies = zeros(length(N),length(B),2);
for i = 1:length(N)
    for j = 1:length(B)
        n=N(i);
        b=B(j);
        strokes = 1:n;

        LPCBX = PCBX(strokes,:);
        LPCBA = PCBA(strokes,:);

        maxLPCBA = max(LPCBA, [], 'all');
        features = transpose(LPCBA) / maxLPCBA;
        tree=fitctree(features,Slabels,'MaxNumSplits',b);

        trainClass = predict(tree, features);
        testClass = predict(tree, transpose(LPCBX) / maxLPCBA);

        treeAccuracies(i,j,1) = length(find(Slabels == trainClass))/length(trainClass);
        treeAccuracies(i,j,2) = length(find(Tlabels == testClass))/length(testClass);
    end
end

%% Visualize accuracies
% Code for plotting values on bar graph modified from Bar Graph
% matlab documentation example.
figure(7)
title('Accuracy of Training and Test Fit')


subplot(1,2,1)
b = bar(treeAccuracies(:,:,1));
xlabel('Number of V-Modes Used'); xticklabels(N)
ylabel('Accuracy of Fit')
axis([0.5 length(N)+0.5 0 1])

xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = strrep(string(round(b(1).YData, 2)), '0.', '.');
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = strrep(string(round(b(2).YData, 2)), '0.', '.');
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips3 = b(3).XEndPoints;
ytips3 = b(3).YEndPoints;
labels3 = strrep(string(round(b(3).YData, 2)), '0.', '.');
text(xtips3,ytips3,labels3,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips4 = b(4).XEndPoints;
ytips4 = b(4).YEndPoints;
labels4 = strrep(string(round(b(4).YData, 2)), '0.', '.');
text(xtips4,ytips4,labels4,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

lgd = legend(string(B));
title(lgd, 'Max Branch Points')
title('Decision Tree Training Accuracy')


subplot(1,2,2)
b = bar(treeAccuracies(:,:,2));
xlabel('Number of V-Modes Used'); xticklabels(N)
ylabel('Accuracy of Fit')
axis([0.5 length(N)+0.5 0 1])

xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = strrep(string(round(b(1).YData, 2)), '0.', '.');
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = strrep(string(round(b(2).YData, 2)), '0.', '.');
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips3 = b(3).XEndPoints;
ytips3 = b(3).YEndPoints;
labels3 = strrep(string(round(b(3).YData, 2)), '0.', '.');
text(xtips3,ytips3,labels3,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips4 = b(4).XEndPoints;
ytips4 = b(4).YEndPoints;
labels4 = strrep(string(round(b(4).YData, 2)), '0.', '.');
text(xtips4,ytips4,labels4,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

lgd = legend(string(B));
title(lgd, 'Max Branch Points')
title('Decision Tree Testing Accuracy')

%% SVM classifier with all digits
PCBX = U'*X; % U'*A = S*V'
PCBA = U'*A;
n = 25;
strokes = 1:n;
% Projections onto V-modes specified by strokes
LPCBX = PCBX(strokes,:);
LPCBA = PCBA(strokes,:);

maxLPCBA = max(LPCBA, [], 'all');
features = transpose(LPCBA) / maxLPCBA;
Mdl = fitcecoc(features,Slabels);

trainClass = predict(Mdl, features);
testClass = predict(Mdl, transpose(LPCBX) / maxLPCBA);

trainAccuracy = length(find(Slabels == trainClass))/length(trainClass)
testAccuracy = length(find(Tlabels == testClass))/length(testClass)

%% Decision Tree Analysis with easy and difficult digits
SdigitMask = find((Slabels == 1) | (Slabels == 4));
TdigitMask = find((Tlabels == 1) | (Tlabels == 4));

n = 25;
strokes = 1:n;

PCBX = U'*X; % U'*A = S*V'
LPCBX = PCBX(strokes,TdigitMask); % Projections onto V-modes specified by strokes
PCBA = U'*A;
LPCBA = PCBA(strokes,SdigitMask);

maxLPCBA = max(LPCBA, [], 'all');
features = transpose(LPCBA) / maxLPCBA;
tree=fitctree(features,Slabels(SdigitMask),'MaxNumSplits',1000);

trainClass = predict(tree, features);
testClass = predict(tree, transpose(LPCBX) / maxLPCBA);

trainAccuracy = length(find(Slabels(SdigitMask) == trainClass))/length(trainClass)
testAccuracy = length(find(Tlabels(TdigitMask) == testClass))/length(testClass)

%% SVM classifier with easy and difficult digits
SdigitMask = find((Slabels == 5) | (Slabels == 8));
TdigitMask = find((Tlabels == 5) | (Tlabels == 8));

PCBX = U'*X; % U'*A = S*V'
PCBA = U'*A;
n = 25;
strokes = 1:n;
% Projections onto V-modes specified by strokes
LPCBX = PCBX(strokes,TdigitMask);
LPCBA = PCBA(strokes,SdigitMask);

maxLPCBA = max(LPCBA, [], 'all');
features = transpose(LPCBA) / maxLPCBA;
Mdl = fitcsvm(features,Slabels(SdigitMask));

trainClass = predict(Mdl, features);
testClass = predict(Mdl, transpose(LPCBX) / maxLPCBA);

trainAccuracy = length(find(Slabels(SdigitMask) == trainClass))/length(trainClass)
testAccuracy = length(find(Tlabels(TdigitMask) == testClass))/length(testClass)



