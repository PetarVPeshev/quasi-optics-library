function af = array_factor(nx, ny, dx, dy, k, k_comp, theta_p, phi_p)
%ARRAY_FACTOR Summary of this function goes here
%   Detailed explanation goes here
    nx = permute(0 : 1 : nx - 1, [1, 3, 2]) ...
        .* ones( [size(k_comp, 1, 2)] );
    ny = permute(0 : 1 : ny - 1, [1, 3, 2]) ...
        .* ones( [size(k_comp, 1, 2)] );

    af_x = exp(1j * dx * nx ...
        .* ( k_comp(:, :, 1) - k * sin(theta_p) .* cos(phi_p) ));
    af_y = exp(1j * dy * ny ...
        .* ( k_comp(:, :, 2) - k * sin(theta_p) .* cos(phi_p) ));
    af_x = permute(sum( permute(af_x, [3, 4, 1, 2]) ), [3, 4, 1, 2]);
    af_y = permute(sum( permute(af_y, [3, 4, 1, 2]) ), [3, 4, 1, 2]);
    
    af = NaN( [size(k_comp, 1, 2), 2] );
    af(:, :, 1) = af_x;
    af(:, :, 2) = af_y;
end

