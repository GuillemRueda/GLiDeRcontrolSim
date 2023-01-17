function F_now = getForce(x_0, model, kc, max_it, beta)
    x_now = x_0 - [0; model.L; 0; 0; 0; 0];
    
    model.rhs = model.d_s - model.M*x_now;
    model.obj = (x_now.')*model.H + model.Alpha;

    V = repmat([model.V1; model.V2], model.N_p, 1);
    S0 = diag(sparse(repmat([1/model.R1; 1/model.R2], model.N_p, 1)));

    curr_results = zeros(9*model.N_p, 1);
    
    for i = 1:max_it
        model.results = curr_results;

        X = model.F*x_now + model.GP*model.results;
        X = [x_now; X(1:(end - 6))];
        X = reshape(X, 6, []);
        r = X(1:3, :);
        r(2, :) = r(2, :) + model.L;
        d = sqrt(sum(r.^2));

        dd = reshape([1./d; sparse(1, model.N_p)], 1, []);
        dd = dd(1:(end - 1));
        
        S = kc*(S0 + diag(dd, 1) + diag(dd, -1));
        
        Q = S\V;
        qq = Q(1:2:(end - 1)).*Q(2:2:end);

        F_C = kc*reshape(r.*repmat(qq.'./(d.^3), 3, 1), [], 1);
        u_C = (model.T - model.kappa_max) * F_C * ...
            (1/model.mass_tug + 1/model.mass_debris);

        model.lb((6*model.N_p + 1):end) = u_C;
        model.ub((6*model.N_p + 1):end) = u_C;

        curr_results = gurobi(model, struct('OutputFlag', 0)).x;
        if sum((curr_results - model.results).^2) < beta
            break;
        end
    end

    model.results = curr_results;

    u_now = model.results(1:3) - model.results(3*model.N_p + (1:3));
    F_now = model.mass_tug*u_now;
end