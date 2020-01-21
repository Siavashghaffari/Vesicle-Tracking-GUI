function outputData = transcytosisIntensity (inputIm, InputCond, TrackCond)
%
% Function is run on a raw (unprocessed) transcytosis video.  Transcytosis
% sites are identified using an adaptive thresholdof the image sequence, 
% with individual transcytosis events determined by the change in
% the slop of the fluorescence intensity of the identified spot.
%
% Processing is done on a single tiff stack; designed to operate at maximum
% speed using parallelization
%
% Detection of exocytosis is via measuing the slope of the vesicle
% intensity across all timepoints of a track except the end points (number
% of end points set by InputCond.finalSteps, min value is 2). A slope
% steper than (InputCond.slopeCutoff) standard deviations the intensity of
% the particle before this point.
%
% Inputs:
%   -inputIm: input file to analyze
%   -InputCond: controls particle detection & graphing
%           -minSize: minimum size (diameter, in pixels) of a detected event to be
%                     classified as trancytosis. 0 = no minimum
%           -maxSize: maximum size (diameter, in pixels) of a detected event to be
%                     classified as trancytosis.  0 = no maximum
%           -circFilter: degree of circularity of a detected object to be
%                        considered an endocytosis event.  Determined using a
%                        measure of circularity, where 1 = perfect circle; 0 =
%                        infinetly long line (i.e. no fiter) [see note 1]
%           -cutOff: minimum intensity a trancytic vessicle must be over the local
%                    background to be scored as a vessicle
%           -psfFWHM: estimated size of the full-width at half-maximumk (FWHM) of
%           the microscopes PSF.  Calculation is:
% 
%               FWHM = ((0.21*wavelength)/NA)/(pixelSize (in um) *1000/mag)
%               i.e. FWHM of 515 emission, 100X microscope, 1.45NA lens, 16um
%               pixels = [(0.21*515)/1.45]/[(16*1000)/100]
%                      = [108.15/1.45]/[16000/100]
%                      = 74.586/160 = 0.46616 pixels
%           -Shift: maximum distance an exosome can move in one frame
%           -TimeVect: time vector for graphs, formatted [start,interval].  End is
%                      calculated automatically
%           -saveGraph: boolean; 1 = save generated graphs 
%           -saveVid: boolean; 1 = save video as colour tiff stack
%           -savgeDiag: boolean; 1 = save diagnostic images
%   -TrackCond: controls particle tracking --> see track.m for full details
%           -mem = number of time steps (frames) that a particle can be
%                  lost and then re-linked
%           -dim = set at 2 unless you know what you're doing
%           -good = minimum track length to keep (in frames)
%           -quiet = 1 suppresses output to command window
%
%  Track.m from: http://physics.georgetown.edu/matlab/index.html
%% default variables to use if not assigned
if nargin == 1
    InputCond.minSize = 3;
    InputCond.maxSize = 9;
    InputCond.circFilter = 0.2;
    InputCond.cutOff = 1.11; %good values - 1.1-1.13???
    InputCond.saveGraph = 1;
    InputCond.saveVid = 1;
    InputCond.saveDiag = 0;
    InputCond.psfFWHM = 0.5;
    InputCond.Shift = 8;
    InputCond.timeVect = [0 0.15];
    InputCond.scale = 0.151498; %scale, in um (or nm) per pixel
    
    TrackCond.mem = 2; %max skip length for linking
    TrackCond.dim = 2; %leave at 2
    TrackCond.good = 5; %min track length
    TrackCond.quiet = 1; %leave at 1
    TrackCond.shortVal = 15; %value below which a track is considered "short"
    TrackCond.shortCol = [0 0 1]; %colour of short tracks on outVid
    TrackCond.longCol = [0 1 0]; %colour of tracks > .shortVal
    
    %requried for intensity analysis
    InputCond.finalSteps = 2; % number of time point to use to look for release
    InputCond.slopeCutoff = 2.5; % number of SD's the final slope must vary from the resot of the slope to be considered exocytosis
end

% Outputs:
%   -outputData: structure containing:
%       -.inputInfo: structure containing input variables, processed
%       folder, and file names
%       -.masks
%           -.mu: cutoffs used to generate thresholds
%           -.mask: mask generated of image variance
%       -.exocytosisSites: structure containing:
%           -.unfiltered: all possible exocyotsis sites
%           -.filtered: all possible exocytosis sites which pass the filter
%       -.frames: strcuture contaiing:
%           -.fileName: name of file processed
%           -.detectedEvents: x/y positions of detected exocytosis events
%       -.graphData: 
%           .timeScale: scale of time, for bottom of graphs
%           .allEvents: table of data for graphs, showing # events per
%            frame, using unfiltered data
%           .filteredEvents: table of data for graphs, showing # events per
%            frame, using filtered data
%   -outVid: TIFF stack of image overlays showing sites of transcytosis.
%
%=================================Notes====================================
% 1) Circularity is calculated as 4*Area/(pi*MaxDiameter^2)
%==========================================================================
%
%  
%

%% check inputs; prepare to process
if (nargin == 0) || (nargin == 2) || (nargin > 3)
    open ('transcytosis.m');
    error ('Incorrect number of input arguments, must call either 1 or 3 arguments.');
end

if (InputCond.circFilter<0) || (InputCond.circFilter>1)
    error ('Valid values for circularity filter are 0 to 1');
end

if (InputCond.minSize<0) || (InputCond.maxSize<0)
    error ('Minimum/Maximum object sizes must be positive intigers');
end

if InputCond.minSize >= InputCond.maxSize
    error ('MaxSize must be larger than minSize');
end

%prepare file lists and output directory

WorkingDir = pwd;
if exist (inputIm, 'file') ~= 2
    error (['Input file ' char(inputIm) ' does not exist.']);
end

a = max(strfind(inputIm, '.'))-1;
folderPre = [inputIm(1:a) '_'];

[~,dirList] = getFileInFolder3('.', folderPre);

if isempty(dirList)
    analysisDir = [folderPre '001'];
else
    numDir = length(dirList)+1;
    a = sprintf('%03d', numDir);
    analysisDir = [folderPre a];
end

mkdir (char(analysisDir));

% clear a numDir forlderPre

% write input info to output file
outputData.inputInfo = InputCond;
outputData.trackCond = TrackCond;
outputData.inputInfo.inputFolder = WorkingDir;
outputData.inputInfo.analyzedFile = inputIm;
outputData.inputInfo.analysisFolder = analysisDir;

%% Prepare for processing

%load stack
tiffStack = LoadStack (inputIm);
stackMin = double(min(tiffStack(:)));
stackMax = double(max(tiffStack(:)));
stack256 = (double(tiffStack) - stackMin)./(stackMax-stackMin);
stack256 = uint8(stack256.*255);
cd (analysisDir);

%generate x-axis titles
maxT = (size(tiffStack,3)*InputCond.timeVect(2))-(InputCond.timeVect(2)-InputCond.timeVect(1));
axisX = [InputCond.timeVect(1):InputCond.timeVect(2):maxT];

%calculate min/max areas to allow
minArea = pi()*(InputCond.minSize/2)^2;
maxArea = pi()*(InputCond.maxSize/2)^2;

if InputCond.minSize == 0
    minArea = 1;
end

if InputCond.maxSize == 0
    maxArea = numel(tiffStack(:,:,1));
end

%regional threshold, filter and analyze image
for i=1:size(tiffStack,3)
    tmpIm = tiffStack(:,:,i);
    ImGauss = imfilter(tmpIm, fspecial('gaussian', [InputCond.minSize InputCond.minSize], 2*InputCond.psfFWHM),'replicate');
    tmpAve = imfilter(ImGauss, fspecial('average', [(2*InputCond.maxSize) (2*InputCond.maxSize)]), 'replicate');
    tmpIm = double(ImGauss)./double(tmpAve);
    ImThresh = false(size(tmpIm)); 
    ImThresh(tmpIm>InputCond.cutOff)=1; 
    outThresh(:,:,i) = ImThresh;
    
    % filter threshold, quantify puncta
    ImStats(i).Unfiltered = regionprops (ImThresh, 'all');
    % clear tmpIm
    tmpIm = bwlabel(ImThresh);
    tmpThresh = false(size(tmpIm));
    statCount = 1; %counter
    for j=1:length(ImStats(i).Unfiltered)
            ImStats(i).Unfiltered(j).Circulatiry = (4*ImStats(i).Unfiltered(j).Area)/(pi()*(ImStats(i).Unfiltered(j).MajorAxisLength^2));
            if (ImStats(i).Unfiltered(j).Area >= minArea) && (ImStats(i).Unfiltered(j).Area <= maxArea && ImStats(i).Unfiltered(j).Circulatiry >= InputCond.circFilter) 
                ImStats(i).Unfiltered(j).Filtered = 1;
                ImStats(i).Filtered(statCount,1) = ImStats(i).Unfiltered(j);
                statCount = statCount+1;
                tmpThresh(tmpIm==j) = 1;
            else
                ImStats(i).Unfiltered(j).Filtered = 0;
            end
    end
    ImFiltThresh(:,:,i) = tmpThresh;
end

outputData.sites = ImStats;

%generate sample variance image and save
write_tiff_stack ('Pre-filter threshold.tif', outThresh);
write_tiff_stack ('Filtered Threshold.tif', ImFiltThresh);
% clear tmpIm tmpAve tmpGauss tmpThresh ImThresh outThresh;

%% Generate data table & track
nParticles = 0;
for i=1:size(tiffStack,3)
    if ~isempty(outputData.sites(i))
        for j=1:length(outputData.sites(i).Filtered)
            nParticles = nParticles + 1;
            outputData.sTable(nParticles,1:2) = outputData.sites(i).Filtered(j).Centroid;
            outputData.sTable(nParticles,3) = (i-1) * InputCond.timeVect(2);
            outputData.sTable(nParticles,4) = outputData.sites(i).Filtered(j).Area; %used to quantify origonal image intensity
            if InputCond.saveDiag > 0
                tTable(i).sites(j,1:2) = outputData.sites(i).Filtered(j).Centroid;
            end
        end
    end %endif
end

save ('transcytosis_analysis.mat', 'outputData');
outputData.tracks = track(outputData.sTable(:,1:3), InputCond.Shift, TrackCond);

%% Generate Track Information
tmpTrack = outputData.tracks;
tmpTrack(:,3) = (tmpTrack(:,3) ./ InputCond.timeVect(2))+1; %convert time to frames
nTracks = max(tmpTrack(:,4));
for i=1:nTracks
    trackx = tmpTrack(tmpTrack(:,4)==i,1);
    tLengths(i) = length(trackx);
    % clear trackx;
end

%% Determine track particle areas
% ID areas for tracks
for i=1:length(outputData.tracks)
    tmpa = find((ismember(outputData.sTable(:,1:3),outputData.tracks(i,1:3))),1,'first'); 
    outputData.tracks(i,5) = outputData.sTable(tmpa,4);
end
% clear tmpa;

for i=1:max(outputData.tracks(:,4))
    tmpArea(i) = mean(outputData.tracks((outputData.tracks(:,4)==i),5));   
end

%% Analyses
% track length histo, velocity histo, MSD & diffusion calcs, 
for i=1:nTracks
    outputData.trackData(i).trackLength = (tLengths(i)-1)*InputCond.timeVect(2); %time in seconds
    outputData.trackData(i).RawPos = tmpTrack(tmpTrack(:,4)==i,1:3);
    outputData.trackData(i).pos = tmpTrack(tmpTrack(:,4)==i,1:2).*InputCond.scale; %convert to scaled units
    outputData.trackData(i).meanArea = tmpArea(i);
    for j=2:size(outputData.trackData(i).pos,1)
        outputData.trackData(i).disp(j-1) = sqrt(((outputData.trackData(i).pos(j,1)-outputData.trackData(i).pos(j-1,1))^2) + ((outputData.trackData(i).pos(j,2)-outputData.trackData(i).pos(j-1,2))^2));
        outputData.trackData(i).vel(j-1) = outputData.trackData(i).disp(j-1)/InputCond.timeVect(2);
    end
    outputData.trackData(i).totalDisp = sum(outputData.trackData(i).disp);
    outputData.trackData(i).meanDisp = mean(outputData.trackData(i).disp);
    outputData.trackData(i).meanVel = mean(outputData.trackData(i).vel);
    [outputData.trackData(i).DcoEff, outputData.trackData(i).Gamma, outputData.trackData(i).MSD] = MSDcalc(outputData.trackData(i).pos, InputCond.timeVect(2));
end
% clear tmpTrack tmpArea

outputData.populationStats.totalTracks = nTracks;
for i=1:nTracks
    outputData.populationStats.trackLength(i) = outputData.trackData(i).trackLength;
    outputData.populationStats.totalDisp(i) = outputData.trackData(i).totalDisp;
    outputData.populationStats.meanDisp(i) = outputData.trackData(i).meanDisp;
    outputData.populationStats.meanVel(i) = outputData.trackData(i).meanVel;
    outputData.populationStats.DcoEff(i) = outputData.trackData(i).DcoEff;
    outputData.populationStats.Gamma(i) = outputData.trackData(i).Gamma;
end

outputData.populationStats.PopMeanTrackLength = mean(outputData.populationStats.trackLength);
outputData.populationStats.PopMeanTotalDisp = mean(outputData.populationStats.totalDisp);
outputData.populationStats.PopMeanDisp = mean(outputData.populationStats.meanDisp);
outputData.populationStats.PopMeanVel = mean(outputData.populationStats.meanVel);
outputData.populationStats.PopDcoEff = mean(outputData.populationStats.DcoEff);
outputData.populationStats.PopGamma = mean(outputData.populationStats.Gamma);
outputData.populationStats.PopMeanTrackLength(2) = std(outputData.populationStats.trackLength)/sqrt(nTracks);
outputData.populationStats.PopMeanTotalDisp(2) = std(outputData.populationStats.totalDisp)/sqrt(nTracks);
outputData.populationStats.PopMeanDisp(2) = std(outputData.populationStats.meanDisp)/sqrt(nTracks);
outputData.populationStats.PopMeanVel(2) = std(outputData.populationStats.meanVel)/sqrt(nTracks);
outputData.populationStats.PopDcoEff(2) = std(outputData.populationStats.DcoEff)/sqrt(nTracks);
outputData.populationStats.PopGamma(2) = std(outputData.populationStats.Gamma)/sqrt(nTracks);

% intensity analysis
possibleExoSites = 0; %counter for exosites to assess further
for i=1:nTracks
    tStrel = strel ('disk', ceil(sqrt(outputData.trackData(i).meanArea/pi())),0);
    tStrel = uint16(tStrel.getnhood);
    strelSize = int16(max(floor(size(tStrel)./2)));
    maxJ = size(outputData.trackData(i).RawPos,1);
    for j=1:maxJ
        tmpIm = tiffStack(:,:,j);
        ImGauss = imfilter(tmpIm, fspecial('gaussian', [InputCond.minSize InputCond.minSize], 2*InputCond.psfFWHM),'replicate');
        
        %edge adjustments for strel
        xStart = int16(1);
        yStart = int16(1);
        xEnd = int16((2*strelSize)+1);
        yEnd = xEnd;
        
        %range on image to quantify
        tmpPos(1) = int16(floor(outputData.trackData(i).RawPos(j,1))-strelSize);
        tmpPos(2) = int16(floor(outputData.trackData(i).RawPos(j,1))+strelSize);
        tmpPos(3) = int16(floor(outputData.trackData(i).RawPos(j,2))-strelSize);
        tmpPos(4) = int16(floor(outputData.trackData(i).RawPos(j,2))+strelSize);
        
        %correct for edges
        if (tmpPos(1) < 1)
            xStart = 1+sum((tmpPos(1):tmpPos(2))<1);
        	tmpPos(1) = 1;
        end
        if ((size(ImGauss,2)) < tmpPos(2))
            xEnd = xEnd - (tmpPos(2)-size(ImGauss,2)); %number of spaces to not use strel
            tmpPos(2) = size(ImGauss,2); %edge of analyzed area
        end
        if (tmpPos(3) < 1)
            yStart = 1 + sum((tmpPos(3):tmpPos(4))<1);
            tmpPos(3) = 1;
        end
        if  ((size(ImGauss,1)) < tmpPos(4))
            yEnd = yEnd - (tmpPos(4)-size(ImGauss,2));
            tmpPos(4) = size(ImGauss,1);
        end
        
        tmpStrel = tStrel (yStart:yEnd, xStart:xEnd);
        tmpVal = ImGauss(tmpPos(3):tmpPos(4),tmpPos(1):tmpPos(2)).*tmpStrel;
        nMeasures = sum(tmpStrel(:));
        outputData.trackData(i).RawPos(j,4) = sum(tmpVal(:))/nMeasures;
         if j == maxJ
             if (max(outputData.trackData(i).RawPos(:,3)) < size(tiffStack,3))
                outputData.trackData(i).RawPos(j+1,1:2) = outputData.trackData(i).RawPos(j,1:2);
                outputData.trackData(i).RawPos(j+1,3) = outputData.trackData(i).RawPos(j,3)+1;
                %measure intensity on next slide
                tmpIm = tiffStack(:,:,j+1);
                ImGauss = imfilter(tmpIm, fspecial('gaussian', [InputCond.minSize InputCond.minSize], 2*InputCond.psfFWHM),'replicate');
                tmpVal = ImGauss(tmpPos(3):tmpPos(4),tmpPos(1):tmpPos(2)).*tmpStrel;
                outputData.trackData(i).RawPos(j+1,4) = sum(tmpVal(:))/nMeasures;
               
             end
             if (max(outputData.trackData(i).RawPos(:,3)) < (size(tiffStack,3)-1))
             %measure intensity on nexter slide
                tmpIm2 = tiffStack(:,:,j+2);
                ImGauss = imfilter(tmpIm2, fspecial('gaussian', [InputCond.minSize InputCond.minSize], 2*InputCond.psfFWHM),'replicate');
                tmpVal = ImGauss(tmpPos(3):tmpPos(4),tmpPos(1):tmpPos(2)).*tmpStrel;
                outputData.trackData(i).RawPos1= sum(tmpVal(:))/nMeasures;
             end
        end % end if j== maxj
        clear tmpVal tmpPos xStart xEnd yStart yEnd tmpIm nMeasures tmpStrel
    end %end for j
    
    %calc fluoro intensity slopes
    for j=2:(size(outputData.trackData(i).RawPos,1)-1) %calc slopes
        outputData.trackData(i).RawPos(j,5) = (outputData.trackData(i).RawPos(j,4)-outputData.trackData(i).RawPos(j-1,4))/(outputData.trackData(i).RawPos(j,3)-outputData.trackData(i).RawPos(j-1,3));
    end
    outputData.trackData(i).RawPos(end,5) = outputData.trackData(i).RawPos(end,4)-outputData.trackData(i).RawPos(end-1,4);
    outputData.trackData(i).RawPos(1,5) = NaN;
    outputData.trackData(i).bleachSlopeMean = mean(outputData.trackData(i).RawPos(2:end-(InputCond.finalSteps),5));
    outputData.trackData(i).bleachSlopeSD = std(outputData.trackData(i).RawPos(2:end-(InputCond.finalSteps),5));
    outputData.trackData(i).finalSlope = outputData.trackData(i).RawPos(end,5);

    
    if (max(outputData.trackData(i).RawPos(:,3)) == size(tiffStack,3) || outputData.populationStats.Gamma(i)>0.87)
        
        outputData.trackData(i).possibleExo = false;
    else
        outputData.trackData(i).possibleExo = true;
        possibleExoSites = possibleExoSites  + 1;
        outputData.possibleExo(possibleExoSites,1) = i;
        outputData.possibleExo(possibleExoSites,2) = outputData.trackData(i).bleachSlopeMean;
        outputData.possibleExo(possibleExoSites,3) = outputData.trackData(i).bleachSlopeSD;
        outputData.possibleExo(possibleExoSites,4) = outputData.trackData(i).finalSlope;
    end
     
    
    clear tStrel strelSize 
end %end for i

%process slopes for real exosites
actSites = 0; %counter to keep track of actual sites
for i=1:size(outputData.possibleExo,1)
	minDif = outputData.possibleExo(i,2) - (InputCond.slopeCutoff * outputData.possibleExo(i,3));
    if (outputData.possibleExo(i,4) <= minDif &&  outputData.trackData(i).RawPos1 <= outputData.trackData(i).RawPos(end,4))
        actSites = actSites + 1;
        outputData.confirmedExo(actSites,:) = outputData.possibleExo(i,:);
    end
end

%generate data for convenient post-processing graphing
%actSites2=2*round(actSites*((actSites^4/100000)+1));
outputData.graphData.nExo = actSites;
outputData.graphData.fReleased = double(actSites)/double(possibleExoSites);
for i=1:actSites
    outputData.graphData.individualExo(i).coords = outputData.trackData(outputData.confirmedExo(i,1)).RawPos(:,1:2); %plot x/y in pixels
    outputData.graphData.individualExo(i).intensity = outputData.trackData(outputData.confirmedExo(i,1)).RawPos(:,3:4); %plot intensity/frame #
    outputData.graphData.individualExo(i).slope(:,1) = outputData.trackData(outputData.confirmedExo(i,1)).RawPos(:,3); %plot slope/frame #
    outputData.graphData.individualExo(i).slope(:,2) = outputData.trackData(outputData.confirmedExo(i,1)).RawPos(:,5);
end
 
%% cleanup
save ('transcytosis_analysis.mat', 'outputData');
cd (WorkingDir);
% clear lastfit;

% display (' ');
% display ('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&');
% display ('&                     Analysis Complete                      &');
% display ('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&');
% display (' ');

end %end function
