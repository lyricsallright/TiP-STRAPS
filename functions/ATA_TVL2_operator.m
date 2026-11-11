function ATAU = ATA_TVL2_operator(U, S, Phi, Kesi, SPC)
%
% ATAop:  ATAU  =  ST FT ( (F S U) Phi )
%   U: [r*Ny*Nx]
%   S: [1, 1, Nc, Ny, Nx]
%   Phi: [r, r, Nc*Ny*Nx]
%
%
% Optional SPC: [] disables spatial constraints; a nonempty mask enables them.
    if nargin < 5
        SPC = [];
    end
    [~,~,Nc,Ny,Nx] = size(S);
    r = size(Phi,1);
    
    % data fidelity term
    if isempty(SPC)
        tmp = reshape(U, [r,1,1,Ny,Nx]);
    else
        U_zero = reshape(U, [r,Ny*Nx]);
        n_background = 1;  % Number of components retained outside the trajectory
        U_zero(n_background+1:end, logical(reshape(SPC,[Ny*Nx,1]))) = 0;
        tmp = reshape(U_zero, [r,1,1,Ny,Nx]);
    end
    tmp = reshape(...
        fft1c_norm(fft1c_norm(S.*tmp, 4), 5), ...
        [r,1,Nc*Ny*Nx]...
    );
    tmp = pagemtimes(Phi, tmp);
    tmp = ifft1c_norm(ifft1c_norm(...
        reshape(tmp, [r,1,Nc,Ny,Nx]), 4), 5);
    tmp = sum(conj(S).*tmp, 3);  % [r,1,1,Ny,Nx]
    
    % TV regularization
    if isempty(SPC)
        tmp2 = Kesi * reshape(U, [r,Ny*Nx]);
    else
        tmp2 = Kesi * (reshape(U, [r,Ny*Nx]) .* reshape(SPC, [1,Ny*Nx]));
    end
    
    % output
    ATAU = tmp(:) + tmp2(:);
    
end

