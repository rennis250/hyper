function [xyY] = XYZ2xyY_rob(XYZ)
xyY = XYZ./sum(XYZ, 2);
xyY(:, 3) = XYZ(:, 2);
return;
end