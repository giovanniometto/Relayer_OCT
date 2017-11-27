function [octVolume, fundus] = readOCTvolumeMEH(hdf5filename)

data_fundus = hdf5read(hdf5filename,'/fundus');
data_oct = hdf5read(hdf5filename,'/oct');
fundus = imrotate(permute(data_fundus,[2 3 1]),-90);
octVolume = permute(data_oct,[2,1,3]);