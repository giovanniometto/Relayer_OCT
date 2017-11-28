function [ILM, RPE, ISOS, THICKNESS] = processVolumeRELAYER(folderORh5file, machineCode, destinationFolder, verbose)

% Author: Giovanni Ometto
% Work address: C274 Tait Building City, University of London, London, EC1V 0HB (UK)
% email: giovanni.ometto@city.ac.uk
% Website: http://www.city.ac.uk
% Jan 2017; Last revision: Oct 2017

% machineCode
% 1: Heidelberg Engineering "Spectralis" px*3.87 = micrometers 
% 2: Topcon "3D oCT-2000" px*2.3 = micrometers 

if nargin<4
    verbose = 0;
end

if ~ischar(folderORh5file)
    disp('Please provide a folder containing sequentially-named image files.')
    return
end

if isdir(folderORh5file)
    
    imExt = ['tif';'jpg';'bmp';'png'];
    tifFiles = dir(fullfile(imFolder,'*.tif'));
    jpgFiles = dir(fullfile(imFolder,'*.jpg'));
    bmpFiles = dir(fullfile(imFolder,'*.bmp'));
    pngFiles = dir(fullfile(imFolder,'*.png'));
    foundExt = [~isempty(tifFiles) ~isempty(jpgFiles) ~isempty(bmpFiles) ~isempty(pngFiles)];
    
    if  sum(foundExt) == 1 
        octVol = makeOctVolumeFromIm (folderORh5file, imExt(find(foundExt),:));
    else
        disp('Please provide a folder containing sequentially-named image files.')
        return
    end
    
else
    
    [~,~,ext] = fileparts(folderORh5file);
    if ext == '.h5'
        data_oct = hdf5read(folderORh5file,'/oct');
        octVol = permute(data_oct,[2,1,3]);
    else
        disp('The provided file is not .h5')
        return
    end

end


if machineCode == 2     % Topcon "3D oCT-2000" px*2.3 = micrometers 
    factor = 0.6693;    % = 2.59/3.87 = 0.6693
    octVol = resizeVolumeOCT(octVol, factor);
else
    octVol = folderORh5file;
end
    
% get parameters
nscans      = size(octVol,3);
scanWidth   = size(octVol,2);

% initialise marices
ILM         = zeros(nscans,scanWidth);
RPE         = zeros(nscans,scanWidth);
ISOS        = zeros(nscans,scanWidth);

% flags
isosUnresolved  = 0;
checkpassed     = 1;

% set 2 start scans as the distance from central scan
distance = 5;
    
if nscans >= (distance*2+1)
    
    mid = round(nscans/2);
    v1= [mid-distance:mid fliplr(1:mid-(distance+1))];
    v2= [fliplr(mid+1:mid+distance) mid mid+(distance+1):nscans];

    for n=1:min(length(v1),length(v2))
            
            if checkpassed

                if n == 1
                    
                    im1 = octVol(:,:,v1(n));
                    im2 = octVol(:,:,v2(n));
                    [ilm1, rpe1, isos1] = segmentFromScratch(im1);
                    ILM(v1(n),:) = ilm1(2,:);
                    RPE(v1(n),:) = rpe1(2,:);
                    ISOS(v1(n),:) = isos1(2,:);
                    ilm1first=fliplr(ilm1');
                    rpe1first=fliplr(rpe1');
                    isos1first=fliplr(isos1');

                    [ilm2, rpe2, isos2] = segmentFromScratch(im2);
                    ILM(v2(n),:) = ilm2(2,:);
                    RPE(v2(n),:) = rpe2(2,:);
                    ISOS(v2(n),:) = isos2(2,:);
                    ilm2first=fliplr(ilm2');
                    rpe2first=fliplr(rpe2');
                    isos2first=fliplr(isos2');

                else


                    if n==distance+1% do the check
                        im = octVol(:,:,v1(n));

                        [ilm1, rpe1, isos1] = segmentFromGuess(im, ilm1, rpe1, isos1);
                        ILM1 = ilm1(2,:);
                        RPE1 = rpe1(2,:);
                        ISOS1 = isos1(2,:);

                        [ilm2, rpe2, isos2] = segmentFromGuess(im, ilm2, rpe2, isos2);
                        ILM2 = ilm2(2,:);
                        RPE2 = rpe2(2,:);
                        ISOS2 = isos2(2,:);

                        madILM = mean(abs(ILM2-ILM1));
                        madRPE = mean(abs(RPE2-RPE1));
                        madISOS = mean(abs(ISOS2-ISOS1));
                        
                        if any([madILM madRPE]>3)
                            
                            checkpassed = 0;
                            % error('An error occurred in the segmentation, cannot continue.')
                            
                            if verbose
                                fprintf( 1, 'check not passed. \n');
                            end
                            
                            [ILM, RPE, ISOS] = alternativeSegmentation(octVol);
                            
                        else
                            if verbose
                                fprintf( 1, 'check passed. proceeding...\n');
                            end

                            if madISOS < 3, isosUnresolved =1; end

                            ILM(v1(n),:) = (ILM1+ILM2)./2;
                            RPE(v1(n),:) = (RPE1+RPE2)./2;
                            ISOS(v1(n),:) = (ISOS1+ISOS2)./2; 
                        end


                    elseif n==2 || n==distance+2 % use the first segmentation
                        im1 = octVol(:,:,v1(n));
                        im2 = octVol(:,:,v2(n));
                        [ilm1, rpe1, isos1] = segmentFromGuess(im1, ilm1first, rpe1first, isos1first);
                        ILM(v1(n),:) = ilm1(2,:);
                        RPE(v1(n),:) = rpe1(2,:);
                        ISOS(v1(n),:) = isos1(2,:);
                        ilm1=fliplr(ilm1');
                        rpe1=fliplr(rpe1');
                        isos1=fliplr(isos1');

                        [ilm2, rpe2, isos2] = segmentFromGuess(im2, ilm2first, rpe2first, isos2first);
                        ILM(v2(n),:) = ilm2(2,:);
                        RPE(v2(n),:) = rpe2(2,:);
                        ISOS(v2(n),:) = isos2(2,:);
                        ilm2=fliplr(ilm2');
                        rpe2=fliplr(rpe2');
                        isos2=fliplr(isos2');

                    else
                        
                        im1 = octVol(:,:,v1(n));
                        im2 = octVol(:,:,v2(n));
                        
                        [ilm1, rpe1, isos1] = segmentFromGuess(im1, ilm1, rpe1, isos1);
                        ILM(v1(n),:) = ilm1(2,:);
                        RPE(v1(n),:) = rpe1(2,:);
                        ISOS(v1(n),:) = isos1(2,:);
                        ilm1=fliplr(ilm1');
                        rpe1=fliplr(rpe1');
                        isos1=fliplr(isos1');

                        [ilm2, rpe2, isos2] = segmentFromGuess(im2, ilm2, rpe2, isos2);
                        ILM(v2(n),:) = ilm2(2,:);
                        RPE(v2(n),:) = rpe2(2,:);
                        ISOS(v2(n),:) = isos2(2,:);
                        ilm2=fliplr(ilm2');
                        rpe2=fliplr(rpe2');
                        isos2=fliplr(isos2');
                        
                    end

                end

            end % checkpassed-conditional block
          
    end
    
    if ~(length(v1)==length(v2)) && checkpassed % last one
        im2 = octVol(:,:,v2(end));
        [ilm2, rpe2, isos2] = segmentFromGuess(im2, ilm2, rpe2, isos2);
        ILM(v2(end),:) = ilm2(2,:);
        RPE(v2(end),:) = rpe2(2,:);
        ISOS(v2(end),:) = isos2(2,:);

    end
    
    if isosUnresolved, priorIsosDistance = (13*size(im,1)/500); ISOS = RPE-priorIsosDistance; end

    THICKNESS = RPE-ILM;
    
    if machineCode == 2
        THICKNESS = THICKNESS./factor;
        RPE = RPE./factor;
        ILM = ILM./factor;
        ISOS = ISOS./factor;
    end
    
    if machineCode == 2
        THICKNESS = THICKNESS.*2.59; % TOPCON 3D-OCT 2000
    else
        THICKNESS = THICKNESS.*3.87; % Heidelberg Engineering Spectralis
    end
 

else % if not 5 or more scans
    
    [ILM, RPE, ISOS] = alternativeSegmentation(octVol);
    
    THICKNESS = RPE-ILM;
    
    if machineCode == 2
        THICKNESS = THICKNESS./factor;
        RPE = RPE./factor;
        ILM = ILM./factor;
        ISOS = ISOS./factor;
    end
    
    if machineCode == 2
        THICKNESS = THICKNESS.*2.59; % TOPCON 3D-OCT 2000
    else
        THICKNESS = THICKNESS.*3.87; % Heidelberg Engineering Spectralis
    end
    
end

if verbose
    surfThickness(THICKNESS, machineCode)
end

% save results
saveResultsAndImages (octVol, destinationFolder, ILM, RPE, ISOS, THICKNESS);
    
end


function [ilm, rpe, isos] = segmentFromScratch(im)

    octMask = getOCTmask(im);

    ilmGuess = getILMguess(im,octMask);
    ilm = snakeILM(im,ilmGuess);
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    rpeGuess  = getRPEguess(im, ilm, octMask);
    rpe = snakeRPE(im,rpeGuess);
   
    isosGuess = getISOSguess(im, rpe, ilm, octMask);
    isos = snakeISOS(im, isosGuess);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    
end


function [ilm, rpe, isos] = segmentFromGuess(im, ilmGuess, rpeGuess, isosGuess)
    ilm = snakeILM(im, ilmGuess);
    rpe = snakeRPE(im, rpeGuess);
    isos = snakeISOS(im, isosGuess);
end


function [ILM, RPE, ISOS] = alternativeSegmentation(octVolume)

fprintf( 1, 'proceeding with alternative segmentation ... \n');

nscans = size(octVolume,3);
mid = round(nscans/2);
v1= (mid):(nscans);
v2= fliplr(1:(mid));

im1 = octVolume(:,:,v1(1));
[ilm1, rpe1, isos1] = segmentFromScratch(im1);
ILM(v1(1),:) = ilm1(2,:);
RPE(v1(1),:) = rpe1(2,:);
ISOS(v1(1),:) = isos1(2,:);

ilm1=fliplr(ilm1');
rpe1=fliplr(rpe1');
isos1=fliplr(isos1');
ilm2=ilm1;
rpe2=rpe1;
isos2=isos1;

for n = 2 : min(length(v2),length(v1))
    
    
    im1 = octVolume(:,:,v1(n));
    im2 = octVolume(:,:,v2(n));
    
    [ilm1, rpe1, isos1] = segmentFromGuess(im1, ilm1, rpe1, isos1);
    ILM(v1(n),:)    = ilm1(2,:);
    RPE(v1(n),:)    = rpe1(2,:);
    ISOS(v1(n),:)   = isos1(2,:);
    ilm1    = fliplr(ilm1');
    rpe1    = fliplr(rpe1');
    isos1   = fliplr(isos1');

    [ilm2, rpe2, isos2] = segmentFromGuess(im2, ilm2, rpe2, isos2);
    ILM(v2(n),:)    = ilm2(2,:);
    RPE(v2(n),:)    = rpe2(2,:);
    ISOS(v2(n),:)   = isos2(2,:);
    ilm2    = fliplr(ilm2');
    rpe2    = fliplr(rpe2');
    isos2   = fliplr(isos2');

end
                    
if ~(length(v1)==length(v2)) % last one
    im1 = octVolume(:,:,v1(end));
    [ilm1, rpe1, isos1] = segmentFromGuess(im1, ilm1, rpe1, isos1);
    ILM(v1(end),:)  = ilm1(2,:);
    RPE(v1(end),:)  = rpe1(2,:);
    ISOS(v1(end),:) = isos1(2,:);
    
end

end


function saveResultsAndImages (octVolume, folder, ILM, RPE, ISOS, THICKNESS)

csvwrite( fullfile(folder,'ILM.csv'),       ILM);
csvwrite( fullfile(folder,'RPE.csv'),       RPE);
csvwrite( fullfile(folder,'ISOS.csv'),      ISOS);
csvwrite( fullfile(folder,'THICKNESS.csv'), THICKNESS);

% save image files of Scans with the segmentation
for i = 1: size(ILM,1)
    
    fh = figure('Visible','off'); imshow(flattenTrimOct(octVolume(:,:,i)),'border','tight');
    
    hold on;
    plot(1:size(ILM,2), ILM(i,:), 'r', 'linewidth',2);
    plot(1:size(RPE,2), RPE(i,:), 'g', 'linewidth',2);    
    hold off
    
    if      i > 99
        saveas(fh,fullfile(folder, [num2str(i) '.jpg'])); % print(fh,fullfile(folder, num2str(i)),'-dbmp256');
    elseif  i > 9
        saveas(fh,fullfile(folder, ['0' num2str(i) '.jpg'])); % print(fh,fullfile(folder, ['0' num2str(i)]),'-dbmp256');
    else
        saveas(fh,fullfile(folder, ['00' num2str(i) '.jpg'])); % print(fh,fullfile(folder, ['00' num2str(i)]),'-dbmp256');
    end
end

end


function surfThickness(THICKNESS, machineCode)

    figure; 
    if machineCode == 2
        px2micron=2.59; % TOPCON 3D-OCT 2000
    else
        px2micron=3.87; % Heidelberg Engineering Spectralis
    end
    figure; 
    surf(THICKNESS.* px2micron, 'edgecolor','none');
    hold on
    axis tight
    colormap jet
    set(gca, 'CLim', [200, 500]);
    cb = colorbar;
    ylabel(cb, 'thickness (\mum)')
    zlim([0 max(THICKNESS(:))*px2micron]);
end


function imoct = flattenTrimOct (imoct)

if size(imoct,3) == 3, imoct = rgb2gray(imoct); end
if (size(imoct,2) > 2*size(imoct,1) && size(imoct,2) ~= 1024), imoct = imoct(:,size(imoct,1)+1:end); end

end
