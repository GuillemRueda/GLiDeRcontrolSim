function model = getModel(m1, m2, L, alpha, d_min, F_max, T, N_p, n, ...
        kappa_max, R1, R2, V1, V2)
    U_max = kappa_max*F_max/m1; % all(U_max >= 0) is true
    
    Phi = getStateTransitionMatrices(T, N_p, n);
    
    F = sparse(blkdiag(Phi{:}))*repmat(sparse(eye(6)), N_p, 1);
    
    G = spalloc(6*N_p, 6*N_p, 36*(N_p+1)*N_p/2);
    for i = 0:(N_p - 1)
        G((1 + 6*i):end, (1 + 6*i):(6 + 6*i)) = ...
            F((1 + 6*i):end, :); %#ok<SPRIX> 
    end
    
    B = repmat({sparse([zeros(3); eye(3)])}, 1, 1, N_p);
    G = G*blkdiag(B{:});
    
    P = sparse([eye(3*N_p), -eye(3*N_p), eye(3*N_p)]);
    Q = sparse(eye(6*N_p));

    model.GP = G*P;
    
    C_s = repmat({sparse([0, -1, 0, 0, 0, 0])}, 1, 1, N_p);
    C_s = blkdiag(C_s{:});
    model.A = C_s*model.GP;
    model.d_s = repmat(-(d_min - L), N_p, 1);
    model.M = C_s*F;
    
    model.Q = model.GP.'*Q*model.GP;
    
    model.H = 2*(F.')*Q*model.GP;
    model.Alpha = repmat(alpha, 1, 9*N_p);
    
    U_max_bar = repmat(U_max(1:3), N_p, 1); % all(U_max_bar >= 0) is true
    U_min_bar = repmat(-U_max(4:6), N_p, 1); % all(U_min_bar <= 0) is true
    
    model.lb = zeros(9*N_p, 1);
    model.ub = [U_max_bar; -U_min_bar; zeros(3*N_p, 1)];
    
    model.L = L;
    model.mass_tug = m1;
    model.mass_debris = m2;
    model.kappa_max = kappa_max;
    model.T = T;
    model.N_p = N_p;
    model.F = F;
    model.R1 = R1;
    model.R2 = R2;
    model.V1 = V1;
    model.V2 = V2;
end