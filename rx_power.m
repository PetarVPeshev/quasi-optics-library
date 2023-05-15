function [rec_power, efficiency] = rx_power(vgo_field, va_field, ...
    rad_power, pw_amplitude, radius, tx_sph_grid, rx_sph_grid)
%RX_POWER Summary of this function goes here
%   Detailed explanation goes here
    wave_impedance = 376.730313668;

    inc_power = pi * (radius ^ 2) * (abs(pw_amplitude) ^ 2) / (2 * wave_impedance);

    theta_max = max(rx_sph_grid(1, end, 1));
    [~, valid_col] = find(tx_sph_grid(:, :, 1) <= theta_max);
    valid_col = unique(valid_col);
    va_field = va_field(:, valid_col, :);

    theta = rx_sph_grid(:, :, 1);
    phi = rx_sph_grid(:, :, 2);

    dth = theta(1, 2) - theta(1, 1);
    dph = phi(2, 1) - phi(1, 1);
    
    va_field = sph2cart_vector(va_field, tx_sph_grid(:, valid_col, :));
    vgo_field = permute(vgo_field, [1 2 4 5 3]);
    vgo_field_x = vgo_field(:, :, :, :, 1) .* sin(theta) .* cos(phi) + ...
        vgo_field(:, :, :, :, 2) .* cos(theta) .* cos(phi) - ...
        vgo_field(:, :, :, :, 3) .* sin(phi);
    vgo_field_y = vgo_field(:, :, :, :, 1) .* sin(theta) .* sin(phi) + ...
        vgo_field(:, :, :, :, 2) .* cos(theta) .* sin(phi) + ...
        vgo_field(:, :, :, :, 3) .* cos(phi);
    vgo_field_z = vgo_field(:, :, :, :, 1) .* cos(theta) - ...
        vgo_field(:, :, :, :, 2) .* sin(theta);
    
    v_oc_dot = va_field(:, :, 1) .* vgo_field_x + va_field(:, :, 2) ...
        .* vgo_field_y + va_field(:, :, 3) .* vgo_field_z;
    v_oc = (2 / wave_impedance) * sum( sum(v_oc_dot .* sin(theta)) ) * dth * dph;
    v_oc = permute(v_oc, [3 4 1 2]);

    rec_power = (abs(v_oc) .^ 2) / (16 * rad_power);
    efficiency = rec_power / inc_power;
end

%     feed_area = pi * (radius ^ 2);
%     tx_sph_grid = tx_sph_grid(:, valid_col, :);
    
%     dth = rx_sph_grid(1, 2, 1) - rx_sph_grid(1, 1, 1);
%     dph = rx_sph_grid(2, 1, 2) - rx_sph_grid(1, 1, 2);

%     vgo_field_cart = NaN( [size(vgo_field)] );

%     vgo_field = sph2cart_vector(vgo_field, rx_sph_grid);

%     v_oc_dot = dot(va_field, vgo_field, 3);

%     inc_power = (abs(pw_amplitude) ^ 2) * feed_area / (2 * wave_impedance);

