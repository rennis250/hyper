clear all; close all; clc;

% This function will take a RAW file from the hyperspectral camera and
% calibrate/convert it. The binary values from the camera will be converted
% to spectral radiance (W/m^2/str/nm). The two outputs of this function are
% the calibrated and converted hyperspectral image (imgSpec) and a vector of
% the measured wavelengths (imgWlns). The image (imgSpec) is a
% 3-dimensional matrix, just like an RGB image, where the dimensions are:
%
% First: height: 800px (always, given how our camera is configured)
% Second: width: depends on how big your scene was and how much you measure
% Third: wavelength measurements: 400 values per pixel always (given how our camera is configured), instead of 3 (RGB)
%
% Note that to properly use this function, you also need the *.hdr file in
% the same directory, so 'apple.hdr' in this case.
[imgSpec, imgWlns] = calibrateRAW('test.raw');

% This will take a hyperspectral image and convert it to various
% colorspaces. The result is a struct (cspaces) that contains fields with
% the colorspaces still in an image format. It needs the hyperspectral image
% and the vector of measured wavelengths, as well as some other info,
% such as the gamma_exps for the monitor that will be used for displaying
% the RGB images. The gamma_exps input should be a vector of 3 values, so
% one gamma exponent for each phosphor.
%
% It then takes a 3x3 matrix that contains the chromaticity coordinates
% (xyY) of the 3 phosphors of the monitor that will be used for
% displaying the RGB images in the following format:
%
% mon_xyY = [x_R, y_R, Y_R;
%            x_G, y_G, Y_G;
%            x_B, y_B, Y_B];
%
% where R, G, and B are the red, green, and blue phosphors of the monitor
% and x, y, and Y are the respective chromaticity coordinates given by a
% device like the Konica Minolta CS2000-A.
%
% Lastly, it expects an illumxyY vector that contains the xyY chromaticity
% coordinates of the illuminant for the scene.
%
% The values I have provided here are for our SONY OLED monitor (the
% gamma_exps and mon_xyY values) and for the illuminant that we used when
% making the fruits/vegetables database (illumxyY).
% gamma_exps = [2.2546, 2.2473, 2.2121];
% mon_xyY = csvread('OLEDxyY.csv');

% But we can also use the sRGB values.
gamma_exps = [2.2, 2.2, 2.2];
mon_xyY = csvread('sRGBxyY.csv');
illumxyY = [0.3324, 0.3435, 36.97];
cspaces = computeColorspaces(imgSpec, imgWlns, gamma_exps, mon_xyY, illumxyY);

% Take a look at the gamma-corrected RGB image.
figure(1);
imshow(cspaces.gRGB);

% Here are the saveCompressedDAT and readCompressedDAT
% functions. These use the first 10 principal components (as extracted by
% the PCA function of Matlab), which contributed roughly 97% of the
% variance for our fruits/vegetables images, which is more than enough for
% color perception purposes (and probably for any spectral purposes).

% Here is how the saveCompressedDAT function is used. Notice that the file
% extension is MAT. The company that makes the camera uses 'DAT' as the
% file extension, so that is why the functions have these names, but we
% save the result in MATLAB's MAT format.
saveCompressedDAT('test.mat', imgSpec, imgWlns);

% Here is how the readCompressedDAT function is used. 
[imgSpec2, imgWlns2] = readCompressedDAT('test.mat');

% We can also apply colormatching functions to the hyperspectral image
% loaded from the MAT file and then gamma correct it for the SONY OLED
% monitor.
cspaces2 = computeColorspaces(imgSpec2, imgWlns2, gamma_exps, mon_xyY, illumxyY);

% We can look at this RGB image to have some idea that saving/loading a compressed
% MAT image does not corrupt the data.
figure(2);
imshow(cspaces2.gRGB);
