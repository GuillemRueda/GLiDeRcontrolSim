function [x_tug, x_debris, delta_v] = simulate(F_now, C, ...
        x_tug, x_debris, mu, m1, m2, F_max, kappa_max, T, t_step, ...
        kc, R1, R2, V1, V2)
    kappa_max = round(kappa_max, 3);
    kappa = round([max(F_now, 0); max(-F_now, 0)]./F_max, 3);
    kappa = min([kappa.'; repmat(kappa_max, 1, 6)]).';
    
    v = unique([0; kappa; kappa_max]);
    l = v(2:end) - v(1:(end-1));
    v = v(1:(end-1));
    n = length(v);
    delta_v = zeros(6, 1);

    U_max = F_max/m1; % all(U_max >= 0) is true

    r_tug = x_tug(1:3);
    v_tug = x_tug(4:6);
    r_debris = x_debris(1:3);
    v_debris = x_debris(4:6);

    for i = 1:n
        tt = l(i);
        u = diag(kappa > v(i))*U_max;
        delta_v = delta_v + u*tt;
        nn = floor(tt/t_step);
        t_vec = [repmat(t_step, nn, 1); tt - nn*t_step];

        for j = 1:(nn + 1)
            t = t_vec(j);
            a_tug = -(mu/(norm(r_tug)^3))*r_tug + [C, -C]*u;
            r_tug = r_tug + v_tug*t + a_tug*(t^2)/2;
            v_tug = v_tug + a_tug*t;
    
            a_debris = -(mu/(norm(r_debris)^3))*r_debris;
            r_debris = r_debris + v_debris*t + a_debris*(t^2)/2;
            v_debris = v_debris + a_debris*t;
        end
    end

    tt = T - kappa_max;
    if tt > 0
        nn = floor(tt/t_step);
        t_vec = [repmat(t_step, nn, 1); tt - nn*t_step];

        for j = 1:(nn + 1)
            t = t_vec(j);
            L = norm(r_tug - r_debris);
            S = kc*[1/R1, 1/L;
                    1/L,  1/R2];
            Q = S\[V1; V2];

            F_c = kc*Q(1)*Q(2)/(L^2);
            F_c = (F_c/L)*(r_tug - r_debris);

            a_tug = -(mu/(norm(r_tug)^3))*r_tug + F_c/m1;
            r_tug = r_tug + v_tug*t + a_tug*(t^2)/2;
            v_tug = v_tug + a_tug*t;
    
            a_debris = -(mu/(norm(r_debris)^3))*r_debris - F_c/m2;
            r_debris = r_debris + v_debris*t + a_debris*(t^2)/2;
            v_debris = v_debris + a_debris*t;
        end
    end

    x_tug = [r_tug; v_tug];
    x_debris = [r_debris; v_debris];
end