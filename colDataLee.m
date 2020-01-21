function ExoData = colDataLee()

% Function to collect data from transcytosis Intensity analysis. 
% Recursively aggregates data from the "transcytosis_analysis.mat into a
% single structure, saved as ExoData.mat. Data is saved in a structure
% called ExoData:
%
%   -.processedDir: The directory that was processed
%   -.fileList: A list of processed files, listed in the same order as the
%               data is processed and recorded into the following data
%               structures.
%   -.populationData: the table I was providing you previously, where the
%                     first column is # exocytic events, second column is
%                     total number of tracked vesicles and the 3rd column
%                     is the fraction of exocytosed tracks. Each row
%                     corrisponds to a single file, with the corrisponding
%                     file having the same row number in fileList.
%   -.individualFiles: a structure containing individual track data for
%                      each file:
%       -.fileName: the file name of the current file
%       -.nExo: Number of detected exocytic events
%       -.nTrack: Number of tracked vesicles
%       -.fExo: Fraction of total tracks that exocytose
%       -.ExoSite: List of x/y coordinates for exocytosis sites
%       -.DockTimes: Time (in frame) a vesicle is docked before exocytosis
%       -.TotalTime: Time (in frames) from the beginning of the recording
%                    until a vesicle is exocytosed
%       -.individualTracks: A structure containing a detailed breakdown of
%                           each track
%           -.trackNumber: Specific track (from the list of all tracks)
%                          which underwent exocytosis
%           -.frames: List of frames each track appears in (in frames, not
%                     units of time)
%           -.coords: x/y coordinates of a vesicle at each time point in 
%                     its track
%           -.intensity: intensity (in picel brightness) of the vesicle at
%                        each time point in its track
%           -.slope: slope of the tracks intensity at each time point in
%                    its track
%           -.size: size of each veiscle (in pixels^2) at each time point
%                   in its track
%
% 
%

%% Find all valid folders, strip extraneous information
homeDir = pwd;
[~,fList] = getFileInFolder ('.');
for ii=1:length(fList)
    if ~isempty(strfind(char(fList(ii)), 'Processed'))
        delList = ii;
    else
        tName = char(fList(ii));
        fList(ii) = cellstr(tName(3:end));
    end
end
%fList(delList) = [];
fList = fList';

clear delList tName

ExoData.processedDir = homeDir;
ExoData.fileList = fList;

%% Collect stats
for ii=1:length(fList);
    cd (char(fList(ii)));
    if exist('transcytosis_analysis.mat', 'file') == 2
        load ('transcytosis_analysis.mat');
        
        %General File Information
        if isfield(outputData, 'confirmedExo')
            
            % collect population stats
            if isfield (outputData, 'graphData')
                ExoData.populationData(ii,1) = outputData.graphData.nExo;
                ExoData.populationData(ii,2) = length(outputData.possibleExo);
                ExoData.populationData(ii,3) = outputData.graphData.fReleased;        
            else
                ExoData.populationData(ii,1) = NaN;
                ExoData.populationData(ii,2) = NaN;
                ExoData.populationData(ii,3) = NaN;
            end
            
            % collect individual track data
            ExoData.individualFiles(ii).fileName = char(fList(ii));
            ExoData.individualFiles(ii).nExo = size(outputData.confirmedExo,1);
            ExoData.individualFiles(ii).nTrack = size (outputData.trackData,2);
            ExoData.individualFiles(ii).fExo = ExoData.individualFiles(ii).nExo/ExoData.individualFiles(ii).nTrack;
            ExoData.individualFiles(ii).ExoSites = []; %variable for exocytosis sites
            ExoData.individualFiles(ii).DockTimes = []; %variable for track dock times
            ExoData.individualFiles(ii).totalTime = []; %variable for track dock times

            for jj=1: size(outputData.confirmedExo,1)
                tNumber = outputData.confirmedExo(jj,1);
                ExoData.individualFiles(ii).individualTrack(jj).trackNumber = tNumber;
                ExoData.individualFiles(ii).individualTrack(jj).frames = outputData.graphData.individualExo(jj).intensity(:,2);
                ExoData.individualFiles(ii).individualTrack(jj).coords = outputData.graphData.individualExo(jj).coords;
                ExoData.individualFiles(ii).individualTrack(jj).intensity = outputData.graphData.individualExo(jj).intensity(:,2);
                ExoData.individualFiles(ii).individualTrack(jj).slope = outputData.graphData.individualExo(jj).slope(:,2);
                ExoData.individualFiles(ii).individualTrack(jj).size = outputData.tracks((outputData.tracks(:,4)==tNumber),5);

                %add data to population stats
                ExoData.individualFiles(ii).ExoSites(jj,1:2) = outputData.graphData.individualExo(jj).coords(end,1:2);
                ExoData.individualFiles(ii).DockTimes(jj) = ExoData.individualFiles(ii).individualTrack(jj).frames(end) - ExoData.individualFiles(ii).individualTrack(jj).frames(1);
                ExoData.individualFiles(ii).totalTime(jj) = ExoData.individualFiles(ii).individualTrack(jj).frames(end);
            end % end jj
        
        else %no exo sites
            ExoData.individualFiles(ii).fileName = char(fList(ii));    
            ExoData.individualFiles(ii).nExo = 0;
            ExoData.individualFiles(ii).nTrack = 0;
            ExoData.individualFiles(ii).fExo = 0;
            ExoData.individualFiles(ii).ExoSites = NaN;
            ExoData.individualFiles(ii).DockTimes = NaN;
            ExoData.populationData(ii,1) = NaN;
            ExoData.populationData(ii,2) = NaN;
            ExoData.populationData(ii,3) = NaN;
        end
            
           
    else
        ExoData.individualFiles(ii).fileName = char(fList(ii));    
        ExoData.individualFiles(ii).nExo = NaN;
        ExoData.individualFiles(ii).nTrack = NaN;
        ExoData.individualFiles(ii).fExo = NaN;
        ExoData.individualFiles(ii).ExoSites = NaN;
        ExoData.individualFiles(ii).DockTimes = NaN;
        ExoData.populationData(ii,1) = NaN;
        ExoData.populationData(ii,2) = NaN;
        ExoData.populationData(ii,3) = NaN;
    end
    if exist ('tracks.fig', 'file')
        open 'tracks.fig'
        saveas (gca, 'tracks.tif');
        close 'all' 'hidden'
    end
    
      cd (homeDir);
end

%% Write Data
cd (homeDir);
save ('ExoData', 'ExoData');



end