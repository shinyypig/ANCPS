function [m0, n0] = IDFT_Integer(im1, im2)
    %the IDFT based translation matching algorithm
    %   only return the integer matching results

    [m, n] = size(im1);
    I1 = fft2((im2));
    I2 = fft2((im1));
    I = I1 .* conj(I2);
    I = I ./ abs(I);
    R = abs(fftshift(ifft2(I)));

    [r, m0] = max(R);
    [~, n0] = max(r);
    m0 = m0(n0);
    m0 = (m - mod(m, 2)) / 2 + 1 - m0;
    n0 = (n - mod(n, 2)) / 2 + 1 - n0;
end
