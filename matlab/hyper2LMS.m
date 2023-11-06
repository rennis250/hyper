function [LMSimg] = hyper2LMS(DAT, wlns)
cssData = csvread('linss2_10e_1.csv');

wavelength_css = cssData(:, 1);
lss = cssData(:, 2);
mss = cssData(:, 3);
sss = cssData(:, 4);

if size(DAT, 3) == 400
    css = interp1(wavelength_css, [lss(:), mss(:), sss(:)], wlns(5:364), 'spline');
else
    % we have a spectrum from a compressed DAT file
    css = interp1(wavelength_css, [lss(:), mss(:), sss(:)], wlns, 'spline');
end

css(css < 0) = 0;

if size(DAT, 3) == 400
    stepSize = [diff(wlns(5:364)); 0];
else
    stepSize = [diff(wlns); 0];
end
corrCSS(:, 1) = css(:, 1) .* stepSize;
corrCSS(:, 2) = css(:, 2) .* stepSize;
corrCSS(:, 3) = css(:, 3) .* stepSize;

DATvis = [];
if size(DAT, 3) == 400
    DATvis = DAT(:, :, 5:364);
else
    DATvis = DAT;
end

h = size(DAT, 1);
w = size(DAT, 2);

DATvis = reshape(DATvis, w*h, size(DATvis, 3));
LMSimg = (corrCSS'*DATvis')';
LMSimg = reshape(LMSimg, h, w, 3);
end