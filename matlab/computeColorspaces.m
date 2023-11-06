function [cspaces] = computeColorspaces(DAT, wlns, gamma_exps, mon_xyY, illumxyY)
[rgb, xyz] = colormatchHyper(DAT, wlns, mon_xyY);
grgb = gammaCorrRGB(rgb, gamma_exps);

hei = size(grgb, 1);
wid = size(grgb, 2);

% separate layers of xyz matrix
Xs = xyz(:, :, 1);
Ys = xyz(:, :, 2);
Zs = xyz(:, :, 3);

% compute xyY coordinates using the layers
s = Xs + Ys + Zs;

xs(:, :) = Xs./s;
ys(:, :) = Ys./s;

xyY(:, :, 1) = xs;
xyY(:, :, 2) = ys;
xyY(:, :, 3) = Ys;

cspaces.xyY = xyY;
cspaces.XYZ = xyz;
cspaces.gRGB = real(grgb);
cspaces.RGB = rgb;

cspaces.LMS = hyper2LMS(DAT, wlns);

illumXYZ = xyY2XYZ_rob(illumxyY);

% cspaces.LAB = xyz2lab(xyz, "WhitePoint", illumXYZ);
lab = xyz2labRob(illumXYZ, [Xs(:), Ys(:), Zs(:)])';
cspaces.LAB = cat(3, reshape(lab(1, :), hei, wid), reshape(lab(2, :), hei, wid), reshape(lab(3, :), hei, wid));

illumX = illumXYZ(1);
illumY = illumXYZ(2);
illumZ = illumXYZ(3);

% now we do LUV, following the formulae on wikipedia
luvls = Ys./illumY;
upn = (4.*illumX)./(illumX + 15.*illumY + 3.*illumZ);
vpn = (9.*illumY)./(illumX + 15.*illumY + 3.*illumZ);

luvycubis = find(luvls > (6/29)^3);
luvylinis = find(luvls <= (6/29)^3);

luvls(luvylinis) = (29/3)^3 .* luvls(luvylinis);
luvls(luvycubis) = 116.*luvls(luvycubis).^(1/3) - 16;

ups = (4.*Xs)./(Xs + 15.*Ys + 3.*Zs);
vps = (9.*Ys)./(Xs + 15.*Ys + 3.*Zs);

luvus = 13.*(luvls) .* (ups - upn);
luvvs = 13.*(luvls) .* (vps - vpn);

luv = [luvls(:), luvus(:), luvvs(:)]';
cspaces.LUV = cat(3, reshape(luv(1, :), hei, wid), reshape(luv(2, :), hei, wid), reshape(luv(3, :), hei, wid));

return
end