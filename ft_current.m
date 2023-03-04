function FT = ft_current(k, k_comp, width, length, antenna, orientation)
%FT_CURRENT This function calculates the current's Fourier Transform
%   Detailed explanation goes here
    keq = k;

    if strcmp(antenna, 'dipole')
        if strcmp(orientation, 'x')
            T = sinc( k_comp(:, :, 2) * width / (2 * pi) );
            F = 2 * keq * ( cos(k_comp(:, :, 1) * length / 2) - ...
                cos(keq * length / 2) ) ./ ( (keq^2 - ...
                k_comp(:, :, 1).^2) * sin(keq * length / 2) );
            FT = zeros( [size(k_comp, 1, 2), 3] );
            FT(:, :, 1) = F .* T;
        else
            error('Not implemented');
        end
    else
        error('Not implemented');
    end
end

