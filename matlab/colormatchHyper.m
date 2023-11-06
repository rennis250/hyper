function [RGBimg, XYZimg] = colormatchHyper(DAT, wlns, mon_xyY)
cmfData = csvread('ciexyz31_1.csv');

wavelength_cmf = cmfData(20:422, 1);

x_bar = cmfData(20:422, 2);
y_bar = cmfData(20:422, 3);
z_bar = cmfData(20:422, 4);

if size(DAT, 3) == 400
    cmf = interp1(wavelength_cmf, [x_bar(:), y_bar(:), z_bar(:)], wlns(5:364), 'spline');
else
    % then we have a compressed DAT file
    cmf = interp1(wavelength_cmf, [x_bar(:), y_bar(:), z_bar(:)], wlns, 'spline');
end

cmf(cmf < 0) = 0;

if size(DAT, 3) == 400
    stepSize = [diff(wlns(5:364)); 0];
else
    stepSize = [diff(wlns); 0];
end
corrCMF(:, 1) = 683 .* cmf(:, 1) .* stepSize;
corrCMF(:, 2) = 683 .* cmf(:, 2) .* stepSize;
corrCMF(:, 3) = 683 .* cmf(:, 3) .* stepSize;

DATvis = [];
if size(DAT, 3) == 400
    DATvis = DAT(:, :, 5:364);
else
    DATvis = DAT;
end

h = size(DAT, 1);
w = size(DAT, 2);

% XYZimg = zeros(h, w, 3);
% for y = 1:h
%     for x = 1:w
%         currBand = squeeze(DATvis(y, x, :))';
%         XYZimg(y, x, :) = currBand*corrCMF;
%     end
% end

DATvis = reshape(DATvis, w*h, size(DATvis, 3));
XYZimg = (corrCMF'*DATvis')';
XYZimg = reshape(XYZimg, h, w, 3);

X = XYZimg(:, :, 1);
Y = XYZimg(:, :, 2);
Z = XYZimg(:, :, 3);

XYZ = [X(:), Y(:), Z(:)];

% monXYZ = xyY2XYZ(mon_xyY);
% RGBimg = xyz2rgb(XYZ, "ColorSpace", "linear-rgb", "WhitePoint", sum(monXYZ, 1));
RGBimg = XYZ2RGB_rob(mon_xyY, XYZ')';
RGBimg = cat(3, reshape(RGBimg(:, 1), h, w), reshape(RGBimg(:, 2), h, w), reshape(RGBimg(:, 3), h, w));
end