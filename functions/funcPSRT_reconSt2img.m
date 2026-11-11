function img = funcPSRT_reconSt2img(reconSt)
    U = reconSt.Us;
    V = reconSt.Vt;
    % data shape
    Nx = size(U,1);
    Ny = size(U,2);
    r = size(U,3);
    Nt = size(V,2);
    % get image
    U = reshape(U, [Nx*Ny,r]);
    img = reshape(abs(U * V), [Nx,Ny,Nt]);
end