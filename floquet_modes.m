function [km_comp, k, modes] = floquet_modes(k0, relat_permit, ...
    dx, dy, sph_grid)
%FLOQUET_MODES Summary of this function goes here
%   Detailed explanation goes here

    k = k0 * sqrt(relat_permit);
    kx0 = k * sin(sph_grid(:, :, 1)) .* cos(sph_grid(:, :, 2));
    ky0 = k * sin(sph_grid(:, :, 1)) .* sin(sph_grid(:, :, 2));
    
    modes = [];
    mx = 0;
    while ~isempty(find((k ^ 2 >= (kx0 - 2 * pi * mx / dx) .^ 2 ...
            + ky0 .^ 2) == 1, 1))
        my = 0;
        while ~isempty(find((k ^ 2 >= (kx0 - 2 * pi * mx / dx) .^ 2 ...
                + (ky0 - 2 * pi * my / dy) .^ 2), 1))
            if mx == 0 && my == 0
                add_modes = [mx my];
            elseif mx == 0 && my ~= 0
                add_modes = [mx my; mx -my];
            elseif mx ~= 0 && my == 0
                add_modes = [mx my; -mx my];
            else
                add_modes = [mx my; mx -my; -mx my; -mx -my];
            end
            modes = [modes; add_modes];
            my = my + 1;
        end
        mx = mx + 1;
    end
    
    % 1D -> phi, 2D -> theta, 3D -> modes
    kxm = 2 * pi * modes(:, 1) / dx;
    kxm = permute(repmat(kxm, [1, 1, size(sph_grid, 1, 2)]), [3 4 1 2]);
    kxm = kx0 - kxm;
    
    kym = 2 * pi * modes(:, 2) / dy;
    kym = permute(repmat(kym, [1, 1, size(sph_grid, 1, 2)]), [3 4 1 2]);
    kym = ky0 - kym;
    
    [~, uhs_col] = find(sph_grid(:, :, 1) <= pi);
    [~, dhs_col] = find(sph_grid(:, :, 1) > pi);

    kzm = NaN( size(kxm) );
    % z >= 0
    kzm(:, uhs_col, :) = - 1j * sqrt( - k .^ 2 ...
        + kxm(:, uhs_col, :) .^ 2 + kym(:, uhs_col, :) .^ 2);
    % z < 0
    kzm(:, dhs_col, :) = 1j * sqrt( - k .^ 2 ...
        + kxm(:, dhs_col, :) .^ 2 + kym(:, dhs_col, :) .^ 2);

    km_comp = NaN( [size(kzm, 1, 2, 3), 3] );
    km_comp(:, :, :, 1) = kxm;
    km_comp(:, :, :, 2) = kym;
    km_comp(:, :, :, 3) = kzm;
end

%     [k_vector, k] = wave_vector(relat_permit, k0, sph_grid);

%     kx0 = k_vector(:, :, 1);
%     ky0 = k_vector(:, :, 2);

%     kx = repmat(k_vector(:, :, 1), [1, 1, 2 * my + 1, 2 * mx + 1]);
    
%     Mx = (-mx : 1 : mx);
%     My = (-my : 1 : my);

%     [Mx, My] = meshgrid((mx : -1 : -mx), (my : -1 : -my));

%     for row_idx = 1 : 1 : length(Mx(:, 1))
%         for col_idx = 1 : 1 : length(Mx(1, :))
%             kxm = kx0 - 2 * pi * Mx(row_idx, col_idx) / dx;
%             kym = ky0 - 2 * pi * My(row_idx, col_idx) / dy;
%             kzm = sqrt( (k ^ 2) - (kxm .^ 2) - (kym .^ 2) );
%         end
%     end
%     ky = repmat(ky, [1, 1, 2 * my]);

%     mode_grid = NaN( [size(Mx, 1, 2), 2] );
%     mode_grid(:, :, 1) = Mx;
%     mode_grid(:, :, 2) = My;
