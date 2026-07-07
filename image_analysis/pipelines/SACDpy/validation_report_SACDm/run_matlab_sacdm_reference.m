function run_matlab_sacdm_reference()
repo = fileparts(fileparts(mfilename('fullpath')));
sacdm_root = '/Users/gmgao/GGscripts/SACDm/SACDm';
out_dir = fullfile(repo, 'validation_report_SACDm', 'matlab_reference');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

addpath(genpath(sacdm_root));
addpath(fullfile(repo, 'validation_report_SACDm', 'matlab_patch'), '-begin');

cases = {
    'left', 560, ...
    fullfile(repo, 'tests', 'testdata', 'TXSHA-undiff-SPEN_JF549_H3K27ac_DL650-SACD-FOV-1-left.tif');
    'right', 647, ...
    fullfile(repo, 'tests', 'testdata', 'TXSHA-undiff-SPEN_JF549_H3K27ac_DL650-SACD-FOV-1-right.tif')
};

for i = 1:size(cases, 1)
    side = cases{i, 1};
    wavelength = cases{i, 2};
    raw_path = cases{i, 3};
    fprintf('Running SACDm for %s (%g nm)\n', side, wavelength);
    stack = read_stack(raw_path);
    result = SACDm(stack, ...
        'pixel', 117, ...
        'wavelength', wavelength, ...
        'NA', 1.45, ...
        'mag', 2, ...
        'iter1', 7, ...
        'iter2', 8, ...
        'ACorder', 2, ...
        'subfactor', 0.8);

    result = single(result);
    write_float_tiff(fullfile(out_dir, ['sacdm-' side '.tif']), result);
    save(fullfile(out_dir, ['sacdm-' side '.mat']), 'result', '-v7');
end
end

function stack = read_stack(path)
info = imfinfo(path);
stack = zeros(info(1).Height, info(1).Width, numel(info), 'single');
for k = 1:numel(info)
    stack(:, :, k) = single(imread(path, k));
end
end

function write_float_tiff(path, image)
t = Tiff(path, 'w');
tagstruct.ImageLength = size(image, 1);
tagstruct.ImageWidth = size(image, 2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 32;
tagstruct.SamplesPerPixel = 1;
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Compression = Tiff.Compression.None;
tagstruct.Software = 'MATLAB SACDm reference runner';
t.setTag(tagstruct);
t.write(image);
t.close();
end
