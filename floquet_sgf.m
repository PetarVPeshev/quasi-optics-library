function SGF = floquet_sgf(relat_permit, k, km_comp)
%FLOQUET_SGF Summary of this function goes here
%   Detailed explanation goes here
    wave_impedance = 376.730313668 / sqrt(relat_permit);

    const = - wave_impedance ./ ( 2 * k * km_comp(:, :, :, 3) );

    SGF = NaN( [size(km_comp, 1, 2, 3), 3, 3] );
    for row = 1 : 1 : 3
        for col = 1 : 1 : 3
            if row == col
                SGF(:, :, :, row, col) = const ...
                    .* (k ^ 2 - km_comp(:, :, :, row) .^ 2);
            else
                SGF(:, :, :, row, col) = const ...
                    .* (- km_comp(:, :, :, row) .* km_comp(:, :, :, col));
            end
        end
    end
end

