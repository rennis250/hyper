function [gr] = gammaCorrRGB(r, gamma_exps)
gr = r;

gr(:, :, 1) = gr(:, :, 1).^(1.0/gamma_exps(1));
gr(:, :, 2) = gr(:, :, 2).^(1.0/gamma_exps(2));
gr(:, :, 3) = gr(:, :, 3).^(1.0/gamma_exps(3));
end