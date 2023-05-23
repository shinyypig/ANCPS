function [m0, n0] = ANCPS(im1, im2)
%the ANCPS algorithm
%   m0, n0 are the decimal matching result

    [m, n] = size(im1);

    I1 = fftshift(fft2((im2)));
    I2 = fftshift(fft2((im1)));
    
    % remove the part that does not participate in the calculation 
    gap = round(min([m, n]) / 4);
    I1 = I1(round(m/2)-gap:round(m/2)+gap, round(n/2)-gap:round(n/2)+gap);
    I2 = I2(round(m/2)-gap:round(m/2)+gap, round(n/2)-gap:round(n/2)+gap);
    
    % calculate NCPS
    I = I1 ./ I2;
    I = I ./ abs(I);
    
    % only keep the components within the circle
    [m_, n_] = size(I);
    [cx, cy] = meshgrid(1:n_, 1:m_);
    mask = sqrt((cx - m_/2 - 0.5).^2 + (cy - n_/2 - 0.5).^2) < min([m_, n_]) / 2;
    I = I .* mask;
    
    % calculate ANCPS
    M = conv2(I, I, "same");
    M = M ./ abs(M);
    
    % generate the vectors that are used to estimate the displacements
    [m_, n_] = size(M);
    [cx, cy] = meshgrid(1:n_, 1:m_);
    mask = sqrt((cx - m_/2 - 0.5).^2 + (cy - n_/2 - 0.5).^2) < min([m_, n_]) / 2 - 1;
    maskx = zeros(m_, n_) > 0;
    maskx(2:end, :) = mask(1:end-1, :);
    masky = zeros(m_, n_) > 0;
    masky(:, 2:end) = mask(:, 1:end-1);
    q = M(mask);
    px = M(maskx);
    py = M(masky);
    
    % estimate the displacements
    m0 = angle(TLS(q(:), px(:))) * m / 2 / pi;
    n0 = angle(TLS(q(:), py(:))) * n / 2 / pi;
end
