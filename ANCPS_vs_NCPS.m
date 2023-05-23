%% The following code was written in MATLAB R2020a.
clear; close all;
addpath(genpath(pwd));

%%
im = imread('cameraman.tif');
% crop the image into the size of 255 * 255 pixels
im1 = double(im(1:end-1, 1:end-1)) / 255;
% generate the cyclic shift matrices
[m, n] = size(im1);
Sx = eye(m);
Sx = Sx(:, [2:end, 1]);
Sy = eye(n);
Sy = Sy([2:end, 1], :);
% cyclically shift the image with (5.5, 5.5) pixels 
im2 = real(Sx^5.5 * im1 * Sy^5.5);

count = 1;
for noise = [0 0.1]
    % add Gaussian noise to the  images
    im1_ = im1 + randn(size(im1)) * noise;
    im2_ = im2 + randn(size(im2)) * noise;
    
    [m, n] = size(im1);
    
    % calculate NCPS
    I1 = fftshift(fft2((im2_)));
    I2 = fftshift(fft2((im1_)));    
    I = I1 ./ I2;
    I = I ./ abs(I);
    
    % calculate ANCPS
    M = conv2(I, I, "same") ./ conv2(ones(m, n), ones(m, n), 'same');
    
    % plot the phase of NCPS and ANCPS
    figure;
    subplot(221), imshow(im1_);
    title('the reference image');
    subplot(222), imshow(im2_);
    title('the cyclically shifted image');
    
    ax3 = subplot(223);
    imshow(angle(I) / pi, []);
    colormap(ax3, parula);
    title('phase of NCPS');
    
    ax4 = subplot(224);
    imshow(angle(M) / pi, []);
    colormap(ax4, parula);
    title('phase of ANCPS');
    if noise == 0
        sgt = sgtitle('noise-free');
        sgt.FontSize = 20;
    else
        sgt = sgtitle('added by Gaussian noise');
        sgt.FontSize = 20;
    end
end
