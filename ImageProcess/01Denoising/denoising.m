close all;clc;

%% Default values
maxI = 256;  % max image intensity


%% Other parameters
Iter = 200;
mu = 0.16;
eta = 0.08;
alpha = 0.01;
beta = 1;
theta = 1;


%% Input image
[Im0, map] = imread("data\trin.bmp");clc


%% Process input image
[Ny, Nx, Nc] = size(Im0);
if Nc > 1; Im0 = rgb2gray(Im0);end
Im0 = double(Im0);
Im0 = (Im0 - min(Im0(:))) / (max(Im0(:)) - min(Im0(:)));


figure(1); imagesc(Im0);colormap("gray");title("Original Image");
%% noise
Im1 = imnoise(Im0,'gaussian',0.251,0.00615);  % gaussian
% Im1 = imnoise(Im0,'speckle',0.2);
Im1 = (Im1 - min(Im1(:))) / (max(Im1(:)) - min(Im1(:)));

% Im0 = Im0 / max(Im0(:));   % Intensity range = [0, 1]

figure(2); imagesc(Im1);colormap("gray");title("Noise Image");

v = Im1;


%% segmentation code
VecParameters = [Iter; mu; eta; alpha; beta; theta];
[u] = active_denosing_mex(single(Im0), single(v), single(VecParameters));
% [u] = testfunction(single(Im0), single(VecParameters));

figure(3);imagesc(u);colormap("gray");title("Denoising Image");








