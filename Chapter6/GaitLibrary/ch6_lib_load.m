function lib = ch6_lib_load(p)
%CH6_LIB_LOAD  Read the gait library from disk, with the checks worth making.
%
%   lib = ch6_lib_load(p)
%
% Reads p.lib.file, or Results/ch6_gait_library.mat when that is empty.
%
% The validation here is not paranoia about file format. A library is used by
% INTERPOLATING between neighbouring entries, and two specific ways of being
% wrong produce a plausible alpha rather than an error:
%
%   * L not ascending -- ch6_lib_alpha brackets with a `find(L <= L_des)` and
%     would interpolate a pair that does not straddle the request, giving a
%     gait for some other step length entirely.
%   * alpha of the wrong SHAPE for the current p -- a library built at
%     bez_deg 5 used with bez_deg 3 would still multiply and add, and produce
%     virtual constraints that are not the ones any gait was solved for.
%
% Both are silent downstream, so they are loud here.
%
% Inputs
%   p : parameter struct (uses p.lib.file, p.ny, p.n_ctrl)
%
% Output
%   lib : the struct written by ch6_lib_build
%
% See also CH6_LIB_BUILD, CH6_LIB_ALPHA, CH6_SIMULATE.

file = p.lib.file;
if isempty(file)
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    file = fullfile(root, 'Results', 'ch6_gait_library.mat');
end

if ~exist(file, 'file')
    error('ch6_lib_load:missing', ...
          ['No gait library at "%s". Build one with\n' ...
           '    lib = ch6_lib_build(ch6_params);\n' ...
           'or set p.lib.file to an existing library.'], file);
end

lib = load(file);

for f = {'L', 'alpha'}
    if ~isfield(lib, f{1})
        error('ch6_lib_load:fields', ...
              'Library "%s" has no .%s (fields: %s).', ...
              file, f{1}, strjoin(fieldnames(lib).', ', '));
    end
end

L = lib.L(:).';
if numel(L) > 1 && any(diff(L) <= 0)
    error('ch6_lib_load:order', ...
          ['Library step lengths are not strictly ascending: %s. ' ...
           'ch6_lib_alpha brackets on this order and would interpolate the ' ...
           'wrong pair.'], mat2str(L, 4));
end

sz = size(lib.alpha);
if numel(sz) < 3, sz(3) = 1; end
if sz(3) ~= numel(L)
    error('ch6_lib_load:count', ...
          'Library has %d step lengths but %d alpha slices.', numel(L), sz(3));
end
if sz(1) ~= p.ny || sz(2) ~= p.n_ctrl
    error('ch6_lib_load:shape', ...
          ['Library alpha is %dx%d but the current parameters want %dx%d ' ...
           '(p.ny = %d outputs, p.bez_deg = %d -> %d coefficients). The ' ...
           'library was built for a different virtual-constraint ' ...
           'parametrization.'], ...
          sz(1), sz(2), p.ny, p.n_ctrl, p.ny, p.bez_deg, p.n_ctrl);
end

lib.file = file;

end
