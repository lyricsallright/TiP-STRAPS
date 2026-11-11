function ATy = AT_operator(y, V, S, M)
%
% ATop:  ATy = ST FT MT y VT
%   y: [Nx, Ny, Nc, Nt]
%   V: [r, Nt]
%   S: [Nx, Ny, Nc]
%   M: [Nx, Ny, Nt]
%
%
    [Nx, Ny, Nc, Nt] = size(y);
    r = size(V, 1);
    M = reshape(M, [Nx,Ny,1,Nt]);
    %
    yy = M .* y;
    yy = reshape(yy, [Nx*Ny*Nc,Nt]) * V';
    yy = reshape(yy, [Nx,Ny,Nc,r]);
    yy = ifft2c_norm(yy);
    yy = sum(conj(S).*yy, 3);  % [Nx, Ny, 1, r]
    ATy = squeeze(yy); 
    
end

