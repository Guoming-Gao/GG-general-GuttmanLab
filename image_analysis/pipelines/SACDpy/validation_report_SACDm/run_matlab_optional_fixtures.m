function run_matlab_optional_fixtures()
repo = fileparts(fileparts(mfilename('fullpath')));
sacdm_root = '/Users/gmgao/GGscripts/SACDm/SACDm';
out_dir = fullfile(repo, 'validation_report_SACDm', 'optional_reference');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

addpath(genpath(sacdm_root));
addpath(fullfile(repo, 'validation_report_SACDm', 'matlab_patch'), '-begin');

%% Wavelet background fixture
[x, y, t] = ndgrid(single(1:128), single(1:96), single(1:3));
bg_input = 18 + 3*sin(x/11) + 2*cos(y/9) + 0.4*t;
bg_input = bg_input + 45*exp(-((x-70).^2 + (y-38).^2) / 180);
bg_input = bg_input + 25*exp(-((x-35).^2 + (y-70).^2) / 100);
bg_output = background_estimation(bg_input, [], 4, 'db6', 2);
save(fullfile(out_dir, 'background_fixture.mat'), 'bg_input', 'bg_output', '-v7');

%% Registration fixture
[xr, yr] = ndgrid(single(1:64), single(1:64));
reg_reference = 0.15 + exp(-((xr-25).^2 + (yr-35).^2) / 50) + 0.65*exp(-((xr-43).^2 + (yr-21).^2) / 28);
reg_reference = single(reg_reference ./ max(reg_reference(:)));
reg_known_shift = [1.35, -2.20];
reg_moving = single(ShiftImage(reg_reference, reg_known_shift));
reg_output = single(register(reg_moving, reg_reference));
save(fullfile(out_dir, 'registration_fixture.mat'), 'reg_reference', 'reg_moving', 'reg_known_shift', 'reg_output', '-v7');

%% Sparse Hessian fixture
[xs, ys] = ndgrid(single(1:18), single(1:16));
sparse_input = single(0.05 + exp(-((xs-8).^2 + (ys-7).^2) / 12) + 0.55*exp(-((xs-13).^2 + (ys-11).^2) / 8));
sparse_output = single(SparseHessian_core(sparse_input, 100, 0.1, 1, 6, 0));
save(fullfile(out_dir, 'sparse_fixture.mat'), 'sparse_input', 'sparse_output', '-v7');
end
