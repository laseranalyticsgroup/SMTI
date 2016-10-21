%% Florian Stroehl    fs417@cam.ac.uk
% 2016-10-21
% v1.0
% 
% - create animations of your SMTI data featuring an outline image,
%     the raw data with mask, the localisations with mask, and an event
%     counter
% - create event histograms showing the rate as events/s and a smoothed 
%     event rate timecourse using a moving average of AVG#events/[X]sec
% - create translation density maps featuring an overlay of the outline
%     image and a 2D histogram of the localisation distribution and a 
%     raw localisation map
% - calculate the event rate and store in/append to Rate.mat


%% clean up workspace and load helper functions
clear all
clc
close all
addpath('./functions/')

%% input
% Imaging parameters
cycleTime = 0.2; % time between two frames [s]
photonBudget = 250; % of a single fluorescent protein per frame (average)
adc2Photon = 13.7; % camera ADC counts per photon
onTime = 0.381; % time between maturation and bleaching [s]
pixelSize = 118; % of the camera [nm]

frameOfTreatment = 300; % 0 if no treatment
framesToRemove = 250; % 0 to use all frames for analysis
maxFrames = 900;

% Used for display in movies only
MinADC = 1800;
MaxADC = 2500;

% Analysis and display parameters
windowSize = 5; % for binning of events [frames]
movingAverageSize = 10; % for moving average of the time course [frames]

% flags
saveMovie = false;
saveRateHistogram = true;
saveDensityMap = true;
saveRate = false;

                   
%% load loc data
warning('off','MATLAB:dlmread:ObsoleteSyntax')
[fileName,pathName,~] = ...
    uigetfile('*.txt','Please select localisation file');
resultPathAndFolderName = pathName; 

disp(['Loading ' pathName fileName '...']);

locData = dlmread([pathName fileName],' ',1);
NLocDataTotal = size(locData,1);
disp(['Number of localisations: ' num2str(NLocDataTotal)]);
NFrames = maxFrames - framesToRemove;

%% apply thresholds
photonThreshold = round(photonBudget*adc2Photon*onTime/cycleTime); 

disp(['Removing localisations below ', ... 
    num2str(photonThreshold), ' ADC...']);
locData(locData(:,4) < photonThreshold, :) = [];
disp(['Removing localisations after frame ', num2str(maxFrames)]);
locData((locData(:,3) > maxFrames),:) = [];
disp(['Removing localisations before frame ', num2str(framesToRemove)]);
locData((locData(:,3) < framesToRemove),:) = [];
NLocData = size(locData,1);
disp(['Number of localisations after thresholding: ' ...
    num2str(NLocData) ' (', ...
    num2str((100*NLocData/NLocDataTotal),'%.1f'), ' %)']);

locDataForRate = (locData - framesToRemove)*cycleTime;

%% load blue channel (BC), MASK, and RAW data
disp('----------------------------------');
BCFileName = [fileName(1:(end-8)) 'BC.tif'];
BC = double(imread([pathName BCFileName]));
BC = BC./max(BC(:));

MASKFileName = [fileName(1:(end-8)) 'BC_MASK.tif'];
MASK = double(imread([pathName MASKFileName]));
MASK = im2bw(MASK);
MASK = bwmorph(MASK,'remove');
MASK = double(MASK);

if saveMovie == true
    RAWFileName = [fileName(1:(end-8)) 'SMTI.tif'];
    RAW = zeros([size(BC),maxFrames-framesToRemove]);


    for i_ = 1:NFrames
        RAW(:,:,i_) = ...
            double(imread([pathName RAWFileName],i_+framesToRemove));
    end 

    RAW(RAW < MinADC) = MinADC;
    RAW(RAW > MaxADC) = MaxADC;
end
[imgh, imgw] = size(BC);

%% create movie
if saveMovie
for i = 1:NFrames
    
    figureHandle1 = figure(1);
    
    % BC in green color
    BCGreen(:,:,1) = zeros(size(BC));
    BCGreen(:,:,2) = BC;
    BCGreen(:,:,3) = zeros(size(BC));
    subplot_tight(2,3,1);
    imagesc(BCGreen)
    axis image
    axis off
    
    % Raw data + white outline
    subplot_tight(2,3,2);
    imagesc(RAW(:,:,i) + MASK*MaxADC)
    colormap(gray)
    caxis([MinADC,MaxADC])
    axis image
    axis off
    
    % white Outline + cumulative localizations
    subplot_tight(2,3,3);
    imagesc(MASK*MaxADC)
    colormap(gray)
    caxis([MinADC,MaxADC])
    axis image
    axis off
    
    hold on
    locX = locData(locData(:,3) < i+framesToRemove,1)/pixelSize;
    locY = locData(locData(:,3) < i+framesToRemove,2)/pixelSize;
    
    plot(locX,locY,'m.')
    
    % rate
    subplot_tight(2,3,[4 5 6],0.08)
    range = ((0:windowSize:(maxFrames-framesToRemove)))*cycleTime;
    locDataTMP = locDataForRate(locDataForRate(:,3)<(i*cycleTime),3);
    HistHandle = histogram(locDataTMP,range);
    axis equal
    axis([0 (maxFrames-framesToRemove)*cycleTime 0 10])
    xlabel 'time [s]'
    ylabel 'events'
    set(gca,'FontName','Helvetica')
    set(gca,'fontsize',12)
    
    % save results
    % grab frame
    F = getframe(figureHandle1);
    try
        imwrite(F.cdata,[pathName fileName(1:(end-9)) ...
            '_movie.tif'],'writemode','append')
    catch
        pause(0.1)
        imwrite(F.cdata,[pathName fileName(1:(end-9)) ...
            '_movie.tif'],'writemode','append')
    end
    
    pause(0.01)
end % i 
end % save?

%% create rate histogram
if saveRateHistogram

% prepare figure
figureHandle2 = figure(2);

% prepare basal and treated data
locDataForHist = (locData - frameOfTreatment)*cycleTime;
locDataForHistPostTreatment = locDataForHist;
locDataForHistPostTreatment(locDataForHistPostTreatment(:,3)<0,:) = [];

% define histogram range
range = ((-frameOfTreatment:windowSize:...
    (maxFrames-frameOfTreatment)))*cycleTime;

% plot histograms
histHandleAll = histogram(locDataForHist(:,3),range);
histDataAll = histHandleAll.Values;
hold on
histHandleTreated = histogram(locDataForHistPostTreatment(:,3),range);
legend('basal','treated')
hold off

% set aspect ratio 2:1:1 (y vs x vs z)
daspect([2 1 1])

% set axis
maxCounts = 10;
frameMin = (framesToRemove-frameOfTreatment)*cycleTime;
frameMax = (maxFrames-frameOfTreatment)*cycleTime;
axis([frameMin frameMax 0 maxCounts])

% label axis
xlabel 'time [s]'
ylabel 'events'

% set font
set(gca,'FontName','Helvetica')
set(gca,'fontsize',12)

% set figure size slightly bigger than plot size to accommodate labels
set(figureHandle2,'units','normalized','position', ... 
    [0.3 0.3 0.5 2.5*maxCounts/(frameMax-frameMin)])

% make axis tight to avoid extra white space in plot
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
axWidth = outerpos(3) - ti(1) - ti(3);
axHeight = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom axWidth axHeight];

% set paper size of saved svg
figureHandle2.PaperPositionMode = 'auto';
figPosition = figureHandle2.PaperPosition;
figureHandle2.PaperSize = [figPosition(3) figPosition(4)];

% save as svg and tif
print(figureHandle2,[pathName fileName(1:(end-9)) ...
            '_RateHistogram'],'-dsvg')
histogramFrame = getframe(figureHandle2);
imwrite(histogramFrame.cdata,...
    [resultPathAndFolderName fileName(1:(end-9)) '_RateHistogram.tif'])

close(figureHandle2)

% smooth histogram with moving average for timecourse
movingAverageKernel = ones(1,movingAverageSize)/movingAverageSize;
histDataAllAveraged = conv(histDataAll,movingAverageKernel,'same');

% prepare figure
figureHandle2 = figure(2);

% plot
maxCounts = 5;
plot(range(1:end-1),histDataAllAveraged,'black-','linewidth',2)
axis([frameMin frameMax 0 maxCounts])

% label axis
xlabel 'time [s]'
ylabel 'average # of events per 10s'

% set font
set(gca,'FontName','Helvetica')
set(gca,'fontsize',12)

% set paper size of saved svg
figureHandle2.PaperPositionMode = 'auto';
figPosition = figureHandle2.PaperPosition;
figureHandle2.PaperSize = [figPosition(3) figPosition(4)];

% save as svg and tif
print(figureHandle2,[pathName fileName(1:(end-9)) ...
            '_RateTimeCourse'],'-dsvg')
rateAverageFrame = getframe(figureHandle2);
imwrite(rateAverageFrame.cdata,...
    [resultPathAndFolderName fileName(1:(end-9)) '_RateTimeCourse.tif'])

close(figureHandle2)

end

%% create density map
if saveDensityMap

% sort localisations into bins of 1x1 um2 size
binsX = round(imgw*pixelSize/1000);
binsY = round(imgh*pixelSize/1000);
dataRoundedY = floor(1 + binsY.*locData(:,2)./(pixelSize*imgh));
dataRoundedX = floor(1 + binsX.*locData(:,1)./(pixelSize*imgw));
binnedData = zeros(binsY,binsX);

for i_ = 1:length(dataRoundedX)
     binnedData(dataRoundedY(i_),dataRoundedX(i_)) = ...
         binnedData(dataRoundedY(i_),dataRoundedX(i_)) + 1;
end

timeSpan = (NFrames - framesToRemove)*cycleTime;
binnedDataRescaled = imresize(binnedData,[imgh,imgw],'cubic')./timeSpan;

% prepare figure
figureHandle3 = figure(3);
if imgh/imgw > 1
    set(figureHandle3,'units','normalized','position', ... 
        [0.3 0.3 0.3*imgw/imgh 0.3])
else
    set(figureHandle3,'units','normalized','position', ... 
        [0.3 0.3 0.3 0.3*imgh/imgw])
end

% plot density map
imagesc(binnedDataRescaled)
axis image
axis off

% plot outline overlay
whiteOutline = cat(3, ones(size(binnedDataRescaled)),...
    ones(size(binnedDataRescaled)),ones(size(binnedDataRescaled)));
hold on 
h = imagesc(whiteOutline); 
hold off
set(h, 'AlphaData', MASK)

% adjust colours
colorbar
colourMax = 0.15;
caxis([0 colourMax])
colormap(morgenstemning())

% stretch from [0 to colourMax]
binnedDataRescaledAndStretched = uint8(256*binnedDataRescaled/colourMax);
binnedDataRescaledAndStretched(MASK == 1) = 255;

% set font
set(gca,'FontName','Helvetica')
set(gca,'fontsize',12)

% set paper size of saved pdf
figureHandle3.PaperPositionMode = 'auto';
figPosition = figureHandle3.PaperPosition;
figureHandle3.PaperSize = [figPosition(3) figPosition(4)];

% save as svg and tif
print(figureHandle3,[pathName fileName(1:(end-9)) ...
            '_DensityMap'],'-dsvg')
imwrite(binnedDataRescaledAndStretched, morgenstemning(), ...
    [resultPathAndFolderName fileName(1:(end-9)) '_DensityMap.tif'])
        
close(figureHandle3)

% sort localisations into bins of 1x1 camera pixel size
dataRoundedY = floor(1 + locData(:,2)./(pixelSize));
dataRoundedX = floor(1 + locData(:,1)./(pixelSize));
roundedData = zeros(imgh,imgw);

for i_ = 1:length(dataRoundedX)
     roundedData(dataRoundedY(i_),dataRoundedX(i_)) = 1;
end

figureHandle3 = figure(4);
if imgh/imgw > 1
    set(figureHandle3,'units','normalized','position', ... 
        [0.3 0.3 0.3*imgw/imgh 0.3])
else
    set(figureHandle3,'units','normalized','position', ... 
        [0.3 0.3 0.3 0.3*imgh/imgw])
end

% blur localisation data slighly
kernel = fspecial('gaussian',[5 5],0.3);
kernel = kernel/max(kernel(:));
smoothLocalisationData = conv2(roundedData,kernel,'same');

% plot localisation map
imagesc(smoothLocalisationData)
cValues = 255;
customColorMap = [linspace(0,1,cValues) 1;...
                  zeros(1,cValues) 1;...
                  linspace(0,1,cValues) 1]';
colormap(customColorMap)
axis image
axis off

% plot outline overlay
whiteOutline = cat(3, ones(size(roundedData)),...
    ones(size(roundedData)),ones(size(roundedData)));
hold on 
h = imagesc(whiteOutline); 
hold off
set(h, 'AlphaData', MASK)

% add mask
smoothLocalisationData = smoothLocalisationData/...
    max(smoothLocalisationData(:));
localisationMapWithMask = uint8(254*smoothLocalisationData);
localisationMapWithMask(MASK == 1) = 255;

% save
imwrite(localisationMapWithMask, customColorMap, ...
    [resultPathAndFolderName fileName(1:(end-9)) '_LocalisationMap.tif'])

close(figureHandle3)

end

%% calculate rate
if saveRate
    if(~exist('Rate.mat','file'))
        % create Rate file
        frameMin = (framesToRemove-frameOfTreatment)*cycleTime;
        frameMax = (maxFrames-frameOfTreatment)*cycleTime;
        rate = frameMin:(windowSize*cycleTime):...
            (frameMax - windowSize*cycleTime);

        % prepare rate file
        save Rate.mat rate
    end

% load rate file
load Rate.mat

% calculate the rate
range = (framesToRemove:windowSize:(maxFrames));
[instant_rate,~] = histcounts(locData(:,3),range);

% append the calculated rate and save
rate = [rate; instant_rate];
save Rate.mat rate
end