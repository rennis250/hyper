function saveCompressedDAT(fn, h, wlns)
hei = size(h, 1);
wid = size(h, 2);

spectra = reshape(h, wid*hei, 400);

% the very short wavelengths have no energy
% in our LED box and often end up being
% very noisy, so we throw those away.
% some of the longer wavelengths are not even
% detected by the human visual system, as far
% as we know, so we throw those away, too.

spectra = spectra(:, 20:364);

% do PCA on the 'cropped' spectra
[w, pc] = pca(spectra', 'Centered', false);

% tests so far show that the first
% seven principal components are enough.
% this is even a bit of overkill according
% to strengths of the returned eigenvalues,
% which are contained in 'ev'
wc = w(:, 1:10);
pcc = pc(:, 1:10);

wlns = wlns(20:364);

% save out the compressed and cleaned data
save(fn, 'wc', 'pcc', 'wid', 'hei', 'wlns');

return;
end