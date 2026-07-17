function [coeffs, x_start] = unpack_z(z, p)
    n1 = p.nu * p.n_coeffs;
    coeffs  = reshape(z(1:n1), [p.nu, p.n_coeffs]);
    x_start = z(n1+1:end);
end
