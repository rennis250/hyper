function caldata = readCAL(spectBin, spatBin)
filename = strcat('Radiometric_calibration_1x1_Scanner.cal');

bands = 800;
samples = 1600;

fid = fopen(filename, 'r');
[data, ~] = fread(fid, bands*samples, 'single');
fclose(fid);

% according to the company, these first two values are always wrong
% and can be corrected to this.
data(1) = 1.63779217004776;
data(2) = 1.6733721792697906;

caldata = reshape(data, [samples, bands]);

if size(caldata, 1) == 1600 && size(caldata, 2) == 800 && (spectBin ~= 1 || spatBin ~= 1)
    binnedcaldata = zeros(round(size(caldata, 1)/spatBin), round(size(caldata, 2)/spectBin));

    spatialC = 1;
    for y = 1:spatBin:size(caldata, 1)
        
        spectralC = 1;
        for wln = 1:spectBin:size(caldata, 2)
            binnedcaldata(spatialC, spectralC) = mean(caldata(y:y+spatBin-1, wln:wln+spectBin-1), 'all');
        
            spectralC = spectralC + 1;
        end
        
        spatialC = spatialC + 1;
    end

    caldata = binnedcaldata;
end
end