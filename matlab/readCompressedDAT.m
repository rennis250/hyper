function [h, wlns] = readCompressedDAT(fn)
load(fn, 'wc', 'pcc', 'wid', 'hei', 'wlns');

spectra = pcc * wc';
h = reshape(spectra', hei, wid, size(spectra, 1));

if numel(wlns) == 400
    wlns = wlns(20:364);
end

return;
end
