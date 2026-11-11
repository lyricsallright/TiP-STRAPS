function reconSt = func_PSRT_recon____DtF_forloop(...
    kdataImg, maskImg, sensMap, V, reconOPT, Spc, U0)

% Optional U0: previous spatial coefficients for a warm start.
    if nargin < 6
        Spc = [];
    end
    if nargin < 7
        U0 = [];
    end
    t_start = tic;
    [Nx, Ny, Nc, Nt] = size(kdataImg);
    r = size(V, 1);
    
    %% prepare vectorized reconstruction operators
    %  ATy: [r,1,1,Ny,Nx]
    ATy = AT_operator(kdataImg, V, sensMap, maskImg);    %(yImg, V, S, M)
    ATy = reshape(permute(ATy,[3,2,1]), [r,1,1,Ny,Nx]);
    ATy = ATy(:);
    %  M: [1,Nt,Ny*Nx]
    M = reshape(permute(maskImg,[3,2,1]), [1,Nt,Ny*Nx]);
    %  Phi: [r,r,Nc*Ny*Nx]
    Vtrans = V.';
    Phi = pagemtimes(Vtrans' .* M,  Vtrans);   %merge-ability
    Phi = reshape(Phi, [r,r,1,Ny,Nx]);
    Phi = repmat(Phi, [1,1,Nc,1,1]);
    Phi = reshape(Phi, [r,r,Nc*Ny*Nx]);
    %  S: [1,1,Nc,Ny,Nx]
    S = reshape(permute(sensMap,[3,2,1]), [1,1,Nc,Ny,Nx]);

    %%  reconstruction settings
    TVcoef = reconOPT.TVcoef;
    % ------------ temporal TV regularization
    % Kesi: TV matrix
    Kesi = eye(Nt);
    Kesi((Nt+1):(Nt+1):end) = -1;
    Kesi(1,1) = 0;
    % TV kernel
    Kesi = V * Kesi;
    Kesi = Kesi * Kesi';
    % transpose (because U is [r,1,1,Ny,Nx])
    Kesi = Kesi.';
    % weighting
    Kesi = TVcoef * Kesi;

    %%  reconstruction
    tol = reconOPT.tol;
    maxiter = reconOPT.maxiter;
    % Normal equation operator
    Mat_op = @(U) ATA_TVL2_operator(U, S, Phi, Kesi, Spc);
    % recon
    if isempty(U0)
        [U, flag, relres, iters] = pcg(...
            Mat_op, ATy, tol, maxiter);
    else
        [U, flag, relres, iters] = pcg(...
            Mat_op, ATy, tol, maxiter, [], [], U0);
    end
    time_all = toc(t_start);
    
    %% output result
    U0 = U;
    U = permute(reshape(U,[r,Ny,Nx]), [3,2,1]);
 
    reconSt = struct(...
        'reconOPT', reconOPT, ...
        'Us',U, 'Vt',V, 'U0', U0, ...
        'flag',flag, 'relres',relres, 'iters',iters, ...
        'time_all', time_all ...
    );

end

