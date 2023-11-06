function [lab] = xyz2labRob(illumXYZ, xyz)
Xs = xyz(:, 1);
Ys = xyz(:, 2);
Zs = xyz(:, 3);

illumX = illumXYZ(1);
illumY = illumXYZ(2);
illumZ = illumXYZ(3);

labnxs = Xs./illumX;
labnys = Ys./illumY;
labnzs = Zs./illumZ;

labxcubis = find(labnxs > (6/29)^3); labxlinis = find(labnxs <= (6/29)^3);
labycubis = find(labnys > (6/29)^3); labylinis = find(labnys <= (6/29)^3);
labzcubis = find(labnzs > (6/29)^3); labzlinis = find(labnzs <= (6/29)^3);

% first we do LAB, following the formulae on wikipedia
linmt = (1/3)*(29/6)^2;
linat = 4/29;

labnxs(labxlinis) = linmt.*(labnxs(labxlinis)) + linat;
labnxs(labxcubis) = labnxs(labxcubis).^(1/3);

labnys(labylinis) = linmt.*(labnys(labylinis)) + linat;
labnys(labycubis) = labnys(labycubis).^(1/3);

labnzs(labzlinis) = linmt.*(labnzs(labzlinis)) + linat;
labnzs(labzcubis) = labnzs(labzcubis).^(1/3);

labls = 116.*labnys - 16;
labas = 500.*(labnxs - labnys);
labbs = 200.*(labnys - labnzs);

lab = [labls(:), labas(:), labbs(:)];
end
