function x = fft1c_norm(x, dim)

x = fftshift(fft(ifftshift(x,dim),[],dim),dim) ...
    / sqrt(size(x,dim));

end
