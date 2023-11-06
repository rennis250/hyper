function [XYZ] = xyY2XYZ_rob(xyY)
% expects monxyY to be in the following format:
% [x_R, y_R, Y_R;
%  x_G, y_G, Y_G;
%  x_B, y_B, Y_B]
%
% get the values from a spectrophotometer
x = xyY(:, 1);
y = xyY(:, 2);
Y = xyY(:, 3);

X = (Y./y) .* x;
Z = (Y./y) .* (1 - y - x);

XYZ = [X(:), Y(:), Z(:)]';
end