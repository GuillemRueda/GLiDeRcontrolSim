function Phi = getStateTransitionMatrices(T, N_p, n)
    t = T*(1:N_p);
    delta_anom = n*t;
    s = sin(delta_anom);
    c = cos(delta_anom);
    
    Phi = zeros(6, 6, N_p);
    Phi(1,1,:) = 4 - 3*c;
    Phi(2,1,:) = 6*(s - delta_anom);
    Phi(4,1,:) = 3*n*s;
    Phi(5,1,:) = 6*n*(c - 1);
    Phi(2,2,:) = 1;
    Phi(3,3,:) = c;
    Phi(6,3,:) = -n*s;
    Phi(1,4,:) = s/n;
    Phi(2,4,:) = 2*(c - 1)/n;
    Phi(4,4,:) = c;
    Phi(5,4,:) = -2*s;
    Phi(1,5,:) = 2*(1 - c)/n;
    Phi(2,5,:) = 4*s/n - 3*t;
    Phi(4,5,:) = 2*s;
    Phi(5,5,:) = 4*c - 3;
    Phi(3,6,:) = s/n;
    Phi(6,6,:) = c;
    
    Phi = mat2cell(Phi, 6, 6, ones(1, N_p));
end