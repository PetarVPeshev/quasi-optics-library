function SGF = dyadic_sgf(relat_permit, k, k_comp, field, current)
%DYADIC_SGF This function calculates the Spectrum-Green-Functions dyad
%   Detailed explanation goes here
    wave_impedance = 376.730313668 / sqrt(relat_permit);

%     k_comp(:, :, 3) = - 1j * sqrt( - k.^2 + ...
%         k_comp(:, :, 1).^2 + k_comp(:, :, 2).^2 );
    const = - wave_impedance ./ ( 2 * k * k_comp(:, :, 3) );
    
    SGF = NaN( [size(k_comp, 1, 2), 3, 3] );
    if strcmp(field, 'E')
        if strcmp(current, 'J')
            for row = 1 : 1 : 3
                for col = 1 : 1 : 3
                    if row == col
                        SGF(:, :, row, col) = const .* ...
                            ( k^2 - k_comp(:, :, row).^2 );
                    else
                        SGF(:, :, row, col) = const .* ...
                            ( - k_comp(:, :, row) .* ...
                            k_comp(:, :, col) );
                    end
                end
            end
        else 
            error('Not implemented');
        end
    else
        error('Not implemented');
    end
end