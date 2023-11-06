function [imgSpec, imgWlns] = calibrateRAW(filename, radScalar, radCorr)
if nargin < 1 || nargin > 3
    error('only 1, 2, or 3 arguments please!')
end

if nargin == 1
    radCorr = readCAL(2, 2);
    radScalar = 1000;
end

if nargin == 2
    radCorr = readCAL(2, 2);
end

[imgSpec, ~, ~, imgWlns, darkFrame, tint, spectralBinning, spatialBinning] = readRAW(filename);

% convert from uW/cm^2/str/nm * radScalar (default = 1000) to W/m^2/str/nm
corrFactor = 1.0e-5/(tint/(radScalar/(spectralBinning*spatialBinning)));

darkSpecAvg = mean(imgSpec(darkFrame+2:end, :, :), 1);
imgSpec = imgSpec(1:darkFrame-1, :, :);

radCorr = corrFactor .* radCorr;
radCorr = reshape(radCorr, 1, size(radCorr, 1), size(radCorr, 2));
imgSpec = bsxfun(@times, bsxfun(@minus, imgSpec, darkSpecAvg), radCorr);

clear radCorr;
clear darkSpecAvg;

imgSpec(imgSpec < 0) = 0;

imgSpec = permute(imgSpec, [2, 1, 3]);

% not sure why, but the hyperspectral images all come out upside down
% now?!...
imgSpec = flipud(imgSpec);

% the first pixel is always corrupted. Specim says there
% is more or less nothing we can do about this.
imgSpec(1) = 0; 
end