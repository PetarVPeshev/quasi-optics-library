function Zin = array_impedance(relat_permit, dipole_l, dipole_w, dx, dy, ...
    k0, k_comp)
%ARRAY_IMPEDANCE Summary of this function goes here
%   Detailed explanation goes here
    % dipoles are oriented along x-axis (SGFxx)
    wave_impedance = 376.730313668 / sqrt(relat_permit);
    k_comp = k_comp * sqrt(relat_permit);
    k = k0 * sqrt(relat_permit);

    mx_max = 40;
    my_max = 40;
    Zin = zeros( [size(k_comp, 1, 2)] );
    for mx = -mx_max : 1 : mx_max
        for my = -my_max : 1 : my_max
            kxm = k_comp(:, :, 1) - 2 * pi * mx / dx;
            kym = k_comp(:, :, 2) - 2 * pi * my / dy;
            kzm = - 1j * sqrt(- k ^ 2 + kxm .^ 2 + kym .^ 2);

            I = 2 * k ...
                * (cos(kxm * dipole_l / 2) - cos(k * dipole_l / 2)) ...
                ./ ((k ^ 2 - kxm .^ 2) * sin(k * dipole_l / 2));
            J = sinc(kym * dipole_w / (2 * pi));
            Gxx = - wave_impedance .* (k ^ 2 - kxm .^ 2) ...
                ./ (2 * k * kzm);

            Zin = Zin + Gxx .* (abs(I) .^ 2) .* (abs(J) .^ 2);
        end
    end
    Zin = - Zin / (dx * dy);
end
% 
% Zin = array_impedance(dipole_l, dipole_w, dx, dy, ...
%     k, km_comp, sgf, k0, k_comp)
%     I = ( 2 * k .* (cos(km_comp(:, :, :, 1) * dipole_l / 2) ...
%         - cos(k * dipole_l / 2)) ) ...
%         ./ ( (k ^ 2 - km_comp(:, :, :, 1) .^ 2) ...
%         * sin(k * dipole_l / 2) );
%     J = sinc( km_comp(:, :, :, 2) * dipole_w / (2 * pi) );
% 
%     prod = permute(sgf(:, :, :, 1, 1) .* (abs(I) .^ 2) ...
%         .* (abs(J) .^ 2), [3, 4, 1, 2]);
%     Zin = permute(- sum(prod) / (dx * dy), [3, 4, 1, 2]);

