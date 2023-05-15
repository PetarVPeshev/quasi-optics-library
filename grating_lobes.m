function lobes = grating_lobes(k, dx, dy, modes)
%GRATING_LOBES Summary of this function goes here
%   Detailed explanation goes here
    theta = linspace(0, 2 * pi, 2001);

    max_mode = max(modes, [], 'all');
    rows = (2 * max_mode + 1) ^ 2;

    kxm = ones(rows, length(theta)) .* cos(theta) * k;
    kym = ones(rows, length(theta)) .* sin(theta) * k;

    mx = kron([1; 0; -1; 1; 0; -1; 1; 0; -1], ones(1, length(theta)));
    my = kron([-1; -1; -1; 0; 0; 0; 1; 1; 1], ones(1, length(theta)));

    kxm = kxm + 2 * pi * mx / dx;
    kym = kym + 2 * pi * my / dy;

    lobes = NaN( [size(kxm, 1, 2), 2] );
    lobes(:, :, 1) = kxm;
    lobes(:, :, 2) = kym;
end

