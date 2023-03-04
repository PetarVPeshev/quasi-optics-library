function field_total = total_field(field)
%TOTAL_FIELD This function calculates the total field magnitude from the
%field's components
%   Detailed explanation goes here
    field_total = sqrt( abs(field(:, :, 1)) .^ 2 ...
        + abs(field(:, :, 2)) .^ 2 + abs(field(:, :, 3)) .^ 2);
end

