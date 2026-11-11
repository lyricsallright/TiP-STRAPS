function x = ifft2c_norm(x)

% short version
x = fftshift(ifft2(ifftshift(x))) * sqrt(size(x,1)*size(x,2));

% % classical version
% x = fftshift(ifft(ifftshift(x,1),[],1),1)*sqrt(size(x,1));
% x = fftshift(ifft(ifftshift(x,2),[],2),2)*sqrt(size(x,2));
end