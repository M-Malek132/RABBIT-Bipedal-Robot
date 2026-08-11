function L = ch5_lie_rows(A, B, c, nmax)
%CH5_LIE_ROWS  Lie derivative rows of an affine function along a linear system.
%
%   L = ch5_lie_rows(A, B, c, nmax)
%
% For xdot = A x + B u and a scalar affine function phi(x) = c x + d, every Lie
% derivative is a matrix power and there is nothing to approximate:
%
%       L_f^i phi        = c A^i x                    (d differentiates away)
%       L_g L_f^(i-1) phi = c A^(i-1) B
%
% RELATIVE DEGREE IS MEASURED, NOT ASSUMED. It is the first i >= 1 at which
% c A^(i-1) B is nonzero. That number is the whole reason the springmass is in
% Chapter 5, and hard-coding "6" would mean a change to the masses, or to which
% cart is constrained, silently produced a controller enforcing the wrong
% derivative. It also means the SAME code handles the relative-degree-1
% velocity constraint used to exercise Section 5.1.
%
% Inputs
%   A, B : n x n, n x nu
%   c    : 1 x n  (the affine part of phi; the constant offset is the caller's)
%   nmax : how many rows to build (>= the expected relative degree + 1)
%
% Output
%   L : struct with
%         .rows   (nmax+1) x n    row i+1 is c A^i,   so row 1 is c
%         .gcoef  (nmax+1) x nu   row i   is c A^(i-1) B
%         .rd     relative degree, or NaN if u never appears within nmax
%
% Cached on (A, B, c, nmax): this sits behind a sampled-data control loop and
% the matrices never change during a run.
%
% See also CH5_BARRIER, CH5_IO_LIN.

persistent keys vals
if isempty(keys), keys = {}; vals = {}; end

key = {A, B, c, nmax};
for i = 1:numel(keys)
    if isequaln(keys{i}, key)
        L = vals{i};
        return;
    end
end

n  = size(A, 1);
nu = size(B, 2);

rows  = zeros(nmax+1, n);
gcoef = zeros(nmax+1, nu);

row = c(:).';
for i = 0:nmax
    rows(i+1, :)  = row;
    gcoef(i+1, :) = row * B;
    row = row * A;
end

% gcoef(i,:) = c A^(i-1) B, so the relative degree is the first nonzero row.
tol = 1e-12 * max(1, norm(gcoef, inf));
rd  = NaN;
for i = 1:nmax+1
    if norm(gcoef(i,:), inf) > tol
        rd = i;
        break;
    end
end

L = struct('rows', rows, 'gcoef', gcoef, 'rd', rd);

keys{end+1} = key; %#ok<AGROW>
vals{end+1} = L;   %#ok<AGROW>

end
