function [par_coeff, per_coeff] = transm_coeff(theta_i, theta_t, ...
    med1_relat_permit, med2_relat_permit)
%TRANSMISSION_COEFF Summary of this function goes here
%   Detailed explanation goes here
    med_1_imped = 376.730313668 / sqrt(med1_relat_permit);
    med_2_imped = 376.730313668 / sqrt(med2_relat_permit);

    per_coeff = 2 * med_1_imped * cos(theta_i) ./ ...
        (med_1_imped * cos(theta_i) + med_2_imped * cos(theta_t));
    par_coeff = 2 * med_1_imped * cos(theta_i) ./ ...
        (med_1_imped * cos(theta_t) + med_2_imped * cos(theta_i));
end

