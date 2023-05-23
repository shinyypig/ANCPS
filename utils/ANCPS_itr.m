function offset = ANCPS_itr(im1, im2, num)
%the iterative verison of ANCPS
%   num is the number of iterations
%   offset is a num*2 matrix, which contains the result of each iteration

    % calculate the integer part of the displacments
%     [mi, ni] = IDFT_Integer(im1, im2);
      loc = IDFT_US(im1, im2, 2);
      mi = round(loc(1));
      ni = round(loc(2));
    
    % crop the two images
    if mi > 0
        im1 = im1(mi+1:end, :);
        im2 = im2(1:end-mi, :);
    elseif mi < 0
        im1 = im1(1:end+mi, :);
        im2 = im2(1-mi:end, :);
    end
    
    if ni > 0
        im1 = im1(:, ni+1:end);
        im2 = im2(:, 1:end-ni);
    elseif mi < 0
        im1 = im1(:, 1:end+ni);
        im2 = im2(:, 1-ni:end);
    end
    
    % initial for iteration
    offset = zeros(num, 2);
    [m, n] = size(im1);
    I2 = fftshift(fft2(im2));
    px = ((1:m) - (m+1)/2)' * ones([1, n]) / m;
    py = ones([m, 1]) * ((1:n) - (n+1)/2) / n;

    mf = 0;
    nf = 0;
    % remove the outermost pixels of the reference image
    im1_ = im1(2:end-1, 2:end-1);
    
    for i = 1:num
        % remove the outermost pixels of the image to be matched
        im2_ = im2(2:end-1, 2:end-1);
        
        % calculate the decimal part of the displacments
        [m_, n_] = ANCPS(im1_, im2_);
        
        % add the new estimated displacements to the cummulative
        % displacements
        mf = m_ + mf;
        nf = n_ + nf;
        
        % cyclically shift the second image
        im2 = ifft2(ifftshift(I2 .* exp(-2j*pi*(mf*px + nf*py))));
        
        % record the estimation of each iteration
        offset(i, 1) = mi + mf;
        offset(i, 2) = ni + nf;
    end
end