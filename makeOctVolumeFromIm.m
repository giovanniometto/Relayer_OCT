function octVolume = makeOctVolumeFromIm (imFolder, ext)

Files = dir(fullfile(imFolder,['*.' ext]));  
nfiles = length(Files);

if nfiles>0

    nSuffix = str2double (Files(1).name(end-6:end-4));
    oct = imread(fullfile(imFolder, Files(1).name));
    oct = flattenTrimOct (oct);
    octVolume = cat(3,oct);
    
    if nfiles > 1
        for i = 2: nfiles

            if str2double (Files(i).name(end-6:end-4)) ~= nSuffix + 1
                disp('Image files do not have sequential suffixes. Cannot continue');
                return
            else
                oct = imread(fullfile(imFolder, Files(i).name));
                oct = flattenTrimOct (oct);
                octVolume = cat(3,octVolume,oct);
                nSuffix = nSuffix+1;
            end

        end
    end
    
else
    disp('No tif images found in the specified directory.');
end


end


function oct = flattenTrimOct (oct)

if size(oct,3) == 3, oct = rgb2gray(oct); end
if size(oct,2) > 2*size(oct,1), oct = oct(:,size(oct,1)+1:end); end

end