close all;clc;

%% Default values



%% Other parameters
Iter = 10;
mu = 0.5;
eta = 0.5;
alpha = 0.01;
beta = 10;
tau = 0.5;
theta = 1;
% labels = 1;


%% Input image
[Im0, map] = imread("data\trin.bmp");


%% Process input image
[Ny, Nx, Nc] = size(Im0);
if Nc > 1; Im0 = rgb2gray(Im0);end
Im0 = double(Im0);
Im0 = Im0 / max(Im0(:));   % Intensity range = [0, 1]
Im0 = Im0 / sum(sum(Im0));
% Im0 = 1 + maxI * Im0;     % Intensity range = [1, maxI]
figure(1); imagesc(Im0);colormap("gray");title("Input image"); colorbar;

c = [30 30 128 128];
r = [30 128 128 30];
BW = roipoly(Im0, c, r);  % u 


%% segmentation code
% gray = (1 + 256) / (labels + 1);

VecParameters = [Iter; mu; eta; alpha; beta; theta; tau;];
[u] = active_segmentation_mex(single(Im0), single(BW), single(VecParameters));
% [u] = testfunction(single(Im0), single(VecParameters));

figure(2);imagesc(u);colormap("gray");title("test");colorbar;








