function [hyper, lines, samples,  wavelength, autodarkFrame, tint, spectralBinning, spatialBinning] = readRAW(filename)
fid = fopen([filename(1:end-4) '.hdr'], 'r');

str = fgets(fid); % read next line
while (str ~= -1) & ~contains(str, 'samples')
    str = fgets(fid);
end

samples = sscanf(str, 'samples = %f');
disp(['samples = ', num2str(samples)]);

lines = sscanf(fgets(fid), 'lines   = %f');
disp(['lines = ', num2str(lines)]);

bands = sscanf(fgets(fid), 'bands = %f');
disp(['bands = ', num2str(bands)]);

autodarkFrame = sscanf(fgets(fid), 'autodarkstartline   = %f');
disp(['autodarkstartline = ', num2str(bands)]);

while (str ~= -1) & ~contains(str, 'binning = {')
    str = fgets(fid);
end

% spatialBinning = fscanf(fid, '%f', 1);
bins = sscanf(str(12:end-1), '%f,', [1, 2]);
spectralBinning = bins(1);
spatialBinning = bins(2);

disp(['spectralBinning = ', num2str(spectralBinning)]);

% spatialBinning = fscanf(fid, '%f}\n', 1);
disp(['spatialBinning = ', num2str(spatialBinning)]);

while (str ~= -1) & ~contains(str, 'tint =  ')
    str = fgets(fid);
end

tint = sscanf(str, 'tint =  %f');
disp(['tint = ', num2str(tint)]);

while (str ~= -1) & ~contains(str, 'Wavelength = {')
    str = fgets(fid);
end

wavelength = fscanf(fid, '%f,\n', bands);

fclose(fid);

hyper = multibandread([filename(1:end-4) '.raw'], [lines, samples, bands], 'uint16', 0, 'bil', 'ieee-le');

% the first value is always infinity and should
% just be set to 0 for simplicity
hyper(1) = 0;
end