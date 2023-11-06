% This function overwhelms my computer. I think it just uses too much
% memory at once. It crashes MATLAB sometimes.

function [imgSpecComp, imgWlnsComp, cspaces] = calibConvCompHyper(fn, gamma_exps, mon_xyY, illumxyY)

warning('This function uses a lot of memory and sometimes crashes MATLAB.');

[imgSpec, imgWlns] = calibrateRAW(fn);

fnwoe = fn(1:end-4);

saveCompressedDAT([fnwoe, '.mat'], imgSpec, imgWlns);
[imgSpecComp, imgWlnsComp] = readCompressedDAT([fnwoe, '.mat']);

cspaces = computeColorspaces(imgSpecComp, imgWlnsComp, gamma_exps, mon_xyY, illumxyY);

return;
end