function [distance, cum_delta_v, delta_vv, ...
    Delta_height, Delta_SMA, i] = execute(T, N_p, kappa_max, ...
          alpha, d_min, L, L0, F_max, N_t, R1, R2, V1, V2, m1, m2, ...
          t_step, max_it, beta, mu, n)

    a = nthroot(mu/(n^2), 3);
    
    x_debris = [a; 0; 0; 0; a*n; 0]; % circular orbit
    x_tug = x_debris + [0; L0; 0; 0; 0; 0];
    
    kc = 8.988e9;
    
    model = getModel(m1, m2, L, alpha, d_min, F_max, T, N_p, n, ...
        kappa_max, R1, R2, V1, V2);
    
    xx_tug = zeros(6, N_t + 1);
    xx_debris = zeros(6, N_t + 1);
    delta_vv = zeros(6, N_t);
    
    xx_tug(:, 1) = x_tug;
    xx_debris(:, 1) = x_debris;
    
    for i = 1:N_t
        try
            [F_now, C] = getThrust(model, x_tug, x_debris, kc, ...
                max_it, beta);
        catch
            break;
        end
        
        [new_x_tug, new_x_debris, delta_v] = simulate(F_now, C, ...
            x_tug, x_debris, mu, m1, m2, F_max, kappa_max, T, t_step, ...
            kc, R1, R2, V1, V2);
    
        xx_tug(:, i + 1) = new_x_tug;
        xx_debris(:, i + 1) = new_x_debris;
        delta_vv(:, i) = delta_v;
    
        x_tug = new_x_tug;
        x_debris = new_x_debris;
    end
    
    distance = sqrt(sum((xx_tug(1:3, :) - xx_debris(1:3, :)).^2));
    cum_delta_v = cumsum(sum(delta_vv));
    Delta_height = sqrt(sum(xx_debris(1:3, :).^2)) - a;
    Delta_SMA = 1./(2./sqrt(sum(xx_debris(1:3, :).^2)) - ...
        (sqrt(sum(xx_debris(4:6, :).^2)).^2)/mu) - a;
end