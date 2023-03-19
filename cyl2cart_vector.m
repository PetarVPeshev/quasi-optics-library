function cart_vector = cyl2cart_vector(cyl_vector, coord_grid)
%CYL2CART_VECTOR This function converts cylindrical vector into cartesian
%   Detailed explanation goes here
    if length( size(coord_grid) ) == 3
        phi = coord_grid(:, :, 2);
    elseif length( size(coord_grid) ) == 2
        phi = coord_grid;
    else
        error('Error. Invalid arguments.');
    end
    
    cart_vector = NaN( [size(cyl_vector, 1, 2), 3] );
    cart_vector(:, :, 1) = cyl_vector(:, :, 1) .* cos(phi) ...
        - cyl_vector(:, :, 2) .* sin(phi);
    cart_vector(:, :, 2) = cyl_vector(:, :, 1) .* sin(phi) ...
        + cyl_vector(:, :, 2) .* cos(phi);
    cart_vector(:, :, 3) = cyl_vector(:, :, 3);
end

