function [F_now, C] = getThrust(model, x_tug, x_debris, kc, max_it, beta)
    r = x_tug(1:3);
    v = x_tug(4:6);
    omega = cross(r, v);
    
    y_coord = v/norm(v);
    z_coord = omega/norm(omega);
    x_coord = cross(y_coord, z_coord);
    C = [x_coord, y_coord, z_coord];
    
    x_rel = [C\(x_tug(1:3) - x_debris(1:3));
        C\(x_tug(4:6) - x_debris(4:6))];
    
    F_now = getForce(x_rel, model, kc, max_it, beta);
end