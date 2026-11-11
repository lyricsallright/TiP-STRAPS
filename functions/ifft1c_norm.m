function x = ifft1c_norm(x, dim)

x = fftshift(ifft(ifftshift(x,dim),[],dim),dim) ...
    * sqrt(size(x,dim));

end
