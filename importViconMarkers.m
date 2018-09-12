function markers = importViconMarkers(varargin)
%IMPORTVICONMARKERS Import numeric data from a text file as a matrix.
%   viconMarkerData = importViconMarkers(path2file,markerNames)
%   NB:This function opens and closes the file twice, the first time to 
%   read the header file and the second time to read the marker positions. 
%   In future it should be re-written to open and close the file once.
%
%   Inputs:: must be in pairs
%       'path2file' - followed by user specified fullfile(path,filename)
%       'markerNames' - cell array containing the names of the markers to
%                       look for
%   Outpurs::
%       markers - struct with fields:
%               .Names - a 'containers.Map' java object whose keys are the
%                        marker names, and whose values denote the columns
%                        in the (.Pos) matrix denote the co-ordinates in R3
%                        of the markers (units in mm)
%               .Pos   - a matrix containing the 3D co-ordinates of
%
%   Example usage::
%   [markers]=importViconMarkers('path2file',path2file,'markerNames',markerNames,'bGaps',bGaps)
%% -------- User Specifies Input args --------------------
bGaps = false;
for i=1:2:nargin
    if  strcmp(varargin{i}, 'path2file'),  path2file = varargin{i+1};
    elseif strcmp(varargin{i}, 'markerNames'),markerNames = varargin{i+1}; 
    elseif strcmp(varargin{i}, 'bGaps'),bGaps = varargin{i+1};
    else error('Invalid argument');
    end    
end

% Initialize variables.
delimiter = ','; % file from vicon will be *.csv
headerBegRow = 3; % get indexes
headerEndRow = 3;

%% Format string for each line of text:
% there will be 3*numberOfMarkers+2 columns
% For more information, see the TEXTSCAN documentation.
[numMarkers,~] = size(markerNames);
numColumns = 2+ numMarkers*3;
headerFormatSpec='';
dataFormatSpec = '';
for cols=1:numColumns
    headerFormatSpec = strcat(headerFormatSpec,'%s');
    dataFormatSpec   = strcat(dataFormatSpec,'%f');
end
headerFormatSpec = strcat(headerFormatSpec,'%[^\n\r]');
dataFormatSpec   = strcat(dataFormatSpec,'%[^\n\r]');
%% Open the text file.
fileIDheader = fopen(path2file,'r');

%% Read header
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileIDheader, '%[^\n\r]', headerBegRow(1)-1, 'ReturnOnError', false);
viconDataHeader = textscan(fileIDheader, headerFormatSpec, headerEndRow(1)-headerBegRow(1)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
viconHeader = [viconDataHeader{1:end-1}];

%% Close the text file.
fclose(fileIDheader);

%% Re-read this time getting ONLY the marker positions and frame numbers
dataBegRow = 6;
dataEndRow = inf;
fileIDmarkers = fopen(path2file,'r');

textscan(fileIDmarkers, '%[^\n\r]', dataBegRow(1)-1, 'ReturnOnError', false);
markerDataArray = textscan(fileIDmarkers, dataFormatSpec, dataEndRow(1)-dataBegRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'ReturnOnError', false);
for block=2:length(dataBegRow)
    frewind(fileIDmarkers);
    textscan(fileIDmarkers, '%[^\n\r]', dataBegRow(block)-1, 'ReturnOnError', false);
    dataArrayBlock = textscan(fileIDmarkers, dataFormatSpec, dataEndRow(block)-dataBegRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'ReturnOnError', false);
    for col=1:length(markerDataArray)
        markerDataArray{col} = [markerDataArray{col};dataArrayBlock{col}];
    end
end
%% Close the text file.
fclose(fileIDmarkers);

%% Create output variable
MarkerPositions = [markerDataArray{1:end-1}];

%% Generate key value pair relationship between axis and column index
numIndexesToFind = length(markerNames);
numColumns       = length(viconHeader);

mapMarkerNames = containers.Map;
for i=1:numIndexesToFind
    for j=1:numColumns
        idx = strfind(viconHeader{j},markerNames{i});
        if ~isempty(idx)
            mapMarkerNames(markerNames{i}) = j:j+2;
        end
    end
end

%% Create output variable
markers        = struct;
markers.Pos    = MarkerPositions;
markers.Names  = mapMarkerNames;

if bGaps
    [markers]=identifyGapsAndPad(markers,markerNames);
end
end

function [markers]=identifyGapsAndPad(markers,markerNames)
    % add additional fields to structure for use without performing gap filling
    markers.GapMarkers = mapMarkerNames.keys; % same fields as markers.Names.keys but sort order is different
    %% ----------- .GapsFilled - logical array denoting frames where markers
    % have dropped out
    [NUM_FRAMES,~] = size(markers.Pos);
    NUM_MARKERS = length(markers.GapMarkers);
    markerGaps = false(NUM_FRAMES,NUM_MARKERS); % assume no gaps initially
    for frameNum=1:NUM_FRAMES
        currFrame = markers.Pos(frameNum,:);
        for m=1:NUM_MARKERS                  
            currMarker=mapMarkerNames(markerNames{m});
            if isnan(currFrame(currMarker));   
                markerGaps(frameNum,m) = true; % gap fill required on this marker at this frame
            end
        end
    end
    markers.GapsFilled = markerGaps;  
    markers.GapsPerFrame = sum(markerGaps,2);
    markers.GapFrames = find(markers.GapsPerFrame~=0);

    % Visual check if missing frames
    vFrames = markers.Pos(:,1);
    frameDiff = vFrames(2:end) - vFrames(1:end-1);
    fhand = figure('name','Check difference of subsequent frames: ');
    plot(frameDiff);

    %% ------------------ Handle Missing Frames if they exist
    if ~isempty(find(frameDiff~=1, 1));
        fprintf('Frames missing\n');
        figure('name','Frame Inter-arrival');plot(frameDiff);
        realFrames = [vFrames(1):vFrames(end)]';
        NUM_REAL_FRAMES = length(realFrames);
        missFrameIdx = true(NUM_REAL_FRAMES,1);
        missFrameIdx(vFrames) = false; % missing frames true, present frames false
        % --- Copy gapFillData to a new struct where missing frames now contain NaN
        markersNanPadData            = struct;
        markersNanPadData.Names      = markers.Names;
        markersNanPadData.GapMarkers = markers.GapMarkers;
        %  --- copy .Pos
        oldPos                  = markers.Pos;
        [~,COLS]                = size(oldPos);
        newPos                  = [realFrames,zeros(NUM_REAL_FRAMES,1),NaN(NUM_REAL_FRAMES,COLS-2)];
        newPos(~missFrameIdx,:) = oldPos;
        markersNanPadData.Pos   = newPos;
        %  --- copy .GapsFilled
        oldGapsFilled                  = markers.GapsFilled;
        [~,COLS2]                      = size(oldGapsFilled);
        newGapsFilled                  = NaN(NUM_REAL_FRAMES,COLS2);
        newGapsFilled(~missFrameIdx,:) = oldGapsFilled;
        markersNanPadData.GapsFilled        = newGapsFilled;
        %  --- copy .GapFrames
        markersNanPadData.GapFrames        = markers.GapFrames;
        %  --- copy .GapsPerFrame
        oldGapsPerFrame                  = markers.GapsPerFrame;
        [~,COLS4]                        = size(oldGapsPerFrame);
        newGapsPerFrame                  = NaN(NUM_REAL_FRAMES,COLS4);
        newGapsPerFrame(~missFrameIdx,:) = oldGapsPerFrame;
        markersNanPadData.GapsPerFrame   = newGapsPerFrame;

        markers = markersNanPadData; % replace original gapFillData struct with Nan padded equivalent
    end
end
