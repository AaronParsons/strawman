Rows = 24;
Cols = 24;
freq = 155*10^6;
Separation = 14;
Latitude = -30.7224 * 2*pi/360;
hourAngle = 0 * 2*pi/24;
timeSteps = 3;
c = 299792458;
totalTimeHours = 100;

saveplots = 1;

hex = 40;
%hex = 0;
rotate = 0;

figure(1); clf

% positions = [];
% if hex > 0
%     for row = hex-1 : -1 : -(hex-1)
%         for col = 0:(2*hex-abs(row))-2
%             xPos = ((-(2*hex-abs(row))+2)/2 + col)*sqrt(3)/2*Separation;
%             yPos = row*Separation*sqrt(3)/2;
%             radius = sqrt(xPos^2 + yPos^2);
%             if radius < 50
%                 positions = [positions; xPos yPos];
%             end
%         end
%     end
%     startOutriggers = length(positions);
%     
%     hex2 = 3;
%     for row = hex2-1 : -1 : -(hex2-1)
%         for col = 0:(2*hex2-abs(row))-2
%             yPos = ((-(2*hex2-abs(row))+2)/2 + col)*sqrt(3)*Separation*(hex-1);
%             xPos = row*Separation*sqrt(3)*(hex-1);
%             radius = sqrt(xPos^2 + yPos^2);
%             if radius < 10000 && ~(xPos == 0 && yPos == 0)
%                 positions = [positions; xPos yPos];
%             end
%         end
%     end 
% else
%     for row = -Rows/2+.5:1:Rows/2-.5
%         for col = -Cols/2+.5:1:Cols/2-.5
%             positions = [positions; row*Separation col*Separation];
%         end
%     end
% end

redundantOutriggers = 1;
positions = [];



for row = hex-1 : -1 : -(hex-1)
    for col = 0:(2*hex-abs(row))-2
        xPos = ((-(2*hex-abs(row))+2)/2 + col)*Separation;
        yPos = row*Separation*sqrt(3)/2;
        radius = sqrt(xPos^2 + yPos^2);
        if radius < 164
            positions = [positions; xPos yPos];
        end
    end
end
numberOfCoreAntennas = length(positions)
startOutriggers = length(positions);

outerPackingRatio = 21;


for row = hex-1 : -1 : -(hex-1)
    for col = 0:(2*hex-abs(row))-2
        xPos = ((-(2*hex-abs(row))+2)/2 + col)*Separation*outerPackingRatio;
        yPos = row*Separation*outerPackingRatio*sqrt(3)/2;
        radius = sqrt(xPos^2 + yPos^2);
        if redundantOutriggers
            if radius < 1100 && radius > 164
                positions = [positions; xPos yPos];
            end
        else
            if radius < 1100 && radius > 164 && (row < 0 || (row == 0 && xPos > 0))
                if yPos < xPos/(sqrt(2)/3) && yPos < -xPos/(sqrt(2)/3) + 100
                    positions = [positions; -xPos -yPos];
                else
                    positions = [positions; xPos yPos];
                end
            end
        end
    end
end
numberOfOutriggers = length(positions) - numberOfCoreAntennas

if rotate
    positionsNew(:,1) = positions(:,2);
    positionsNew(:,2) = positions(:,1);
    positions = positionsNew;
end

longestBaseline = max((positions(1:startOutriggers,1).^2 + positions(1:startOutriggers,2).^2).^.5);

figure(1);
plot(positions(:,1),positions(:,2),'ko','MarkerSize',4)
axis equal
xlabel('East-West Position (m)')
ylabel('North-South Position (m)')
plotSize = 1200;
set(gca,'XLim',[-plotSize plotSize],'YLim',[-plotSize plotSize])
set(1,'Color',[1 1 1])
if saveplots
    if hex > 0
        export_fig HexagonalHERA_Positions.pdf  -nocrop
    else
        export_fig SquareHERA_Positions.pdf  -nocrop
    end
end

%%
baselineMap = containers.Map();
for i = 1:length(positions)
    disp(['Now working on ' num2str(i) '...'])
    for j = 1:length(positions)
        if i ~= j && ~(i > startOutriggers && j > startOutriggers)
            deltaEast = positions(i,1) - positions(j,1);
            deltaWest = positions(i,2) - positions(j,2);
            key = [num2str(deltaEast) ',' num2str(deltaWest)];
            if isKey(baselineMap,key)
                baselineMap(key) = baselineMap(key) + 1;
            else
                baselineMap(key) = 1;
            end
        end
    end
end

%%

baselineList = {};
for f = 7:7
    baselineKeys = keys(baselineMap);
    baselineValues = values(baselineMap);
    baselineList{f} = [];
    for k = 1:length(baselineKeys)
        commaPos = strfind(baselineKeys{k},',');
        deltaEast = str2num(baselineKeys{k}(1:(commaPos-1)));
        deltaNorth = str2num(baselineKeys{k}((commaPos+1):end));
        baselineList{f} = [baselineList{f}; deltaEast deltaNorth];
    end
end

figure(555); close(555); hfig = figure(555); clf;
set(555,'Position',[2007         210         867         723])
lambda = c/150000000;
plot(baselineList{7}(:,1)/lambda, baselineList{7}(:,2)/lambda,'.')
title([num2str((80+10*7)) ' MHz'],'FontSize',14)
xlabel('u'); ylabel('v');
changeFontSize(14);
axis square



%%

uLabelList = {};
thetaLabelList = {};
observationTimesList = {};


for f = 1:12
    freq = (80+10*f)*10^6;
    disp(['Now working on ' num2str(freq/10^6) ' MHz...'])
    FWHM = 9.8 * (freq/(150*10^6))^-1; %Degrees
    
    %angRes = .2*c / freq / (longestBaseline);
    angRes = 150000000/freq * 1/1200;
    
    nx = round(FWHM/360*2*pi / angRes / 2)*2;
    FoV = nx*angRes;
    startHA = hourAngle - FoV/2;
    endHA = hourAngle + FoV/2;
    deltaHA = (endHA - startHA)/timeSteps;
    
    
    baselineKeys = keys(baselineMap);
    baselineValues = values(baselineMap);
    
    discardedTime = 0;
    observationTimes = zeros(nx);
    deltaU = 1/angRes/nx;
    for ha = startHA + deltaHA/2 : deltaHA : endHA - deltaHA/2
        for k = 1:length(baselineKeys)
            commaPos = strfind(baselineKeys{k},',');
            deltaEast = str2num(baselineKeys{k}(1:(commaPos-1)));
            deltaNorth = str2num(baselineKeys{k}((commaPos+1):end));
            
            EarthVector = [0 -sin(Latitude) cos(Latitude); 1 0 0; 0 cos(Latitude) sin(Latitude)] * [deltaEast; deltaNorth; 0]; %Delta Up assumed to be 0;
            dec = Latitude;
            uvwVector = [sin(ha) cos(ha) 0; -sin(dec)*cos(ha) sin(dec)*sin(ha) cos(dec); cos(dec)*cos(ha) -cos(dec)*sin(ha) sin(dec)] * EarthVector / (c / freq);        
            
            uAboveRatio = 1-(uBinAbove - (uvwVector(1)/deltaU + nx/2 + 1));
            vAboveRatio = 1-(vBinAbove - (uvwVector(2)/deltaU + nx/2 + 1));
            
            if ((uBinAbove>0) && (uBinAbove<=nx) && (vBinAbove > 0) && (vBinAbove <= nx))
                observationTimes(vBinAbove,uBinAbove) = observationTimes(vBinAbove,uBinAbove) + uAboveRatio*vAboveRatio*totalTimeHours/timeSteps*60*60 *baselineMap(baselineKeys{k});
            end
            if ((uBinAbove>0) && (uBinAbove<=nx) && (vBinBelow > 0) && (vBinBelow <= nx))
                observationTimes(vBinBelow,uBinAbove) = observationTimes(vBinBelow,uBinAbove) + uAboveRatio*(1-vAboveRatio)*totalTimeHours/timeSteps*60*60 *baselineMap(baselineKeys{k});
            end
            if ((uBinBelow>0) && (uBinBelow<=nx) && (vBinAbove > 0) && (vBinAbove <= nx))
                observationTimes(vBinAbove,uBinBelow) = observationTimes(vBinAbove,uBinBelow) + (1-uAboveRatio)*vAboveRatio*totalTimeHours/timeSteps*60*60 *baselineMap(baselineKeys{k});
            end
            if ((uBinBelow>0) && (uBinBelow<=nx) && (vBinBelow > 0) && (vBinBelow <= nx))
                observationTimes(vBinBelow,uBinBelow) = observationTimes(vBinBelow,uBinBelow) + (1-uAboveRatio)*(1-vAboveRatio)*totalTimeHours/timeSteps*60*60 *baselineMap(baselineKeys{k});
            end
            
            %             uBin = round(uvwVector(1)/deltaU) + nx/2 + 1;
            %             vBin = round(uvwVector(2)/deltaU) + nx/2 + 1;
            %
            %             if ((uBin>0) && (uBin<=nx) && (vBin > 0) && (vBin <= nx))
            %                 observationTimes(vBin,uBin) = observationTimes(vBin,uBin) + totalTimeHours/timeSteps*60*60 *baselineMap(baselineKeys{k});
            %             end
            %
        end
    end
    
    uLabels = -nx/2*deltaU:deltaU:(nx/2*deltaU-deltaU);
    uLabelList{f} = uLabels;
    thetaLabelList{f} = -nx/2*angRes/2/pi*360:angRes/2/pi*360:(nx/2*angRes/2/pi*360-angRes/2/pi*360);
    observationTimesList{f} = observationTimes;
end

%%

figure(55); close(55); hfig = figure(55); clf;
set(55,'Position',[2059          89         678         859]);
handle = tight_subplot(4,3,[.07 .07],[.15 .05],[.06 .01]);
for f = 1:12
    axes(handle(f));
    imagesc(uLabelList{f},uLabelList{f},log10(observationTimesList{f}))
    colormap jet
    title([num2str((80+10*f)) ' MHz'],'FontSize',14)
    xlabel('u'); ylabel('v');
    caxis([5 8.5])
end

axes(handle(1));
pos1 = get(handle(1),'Position');
clrbar = colorbar('location','southoutside');
set(handle(1),'Position',pos1);

pos12 = get(handle(12),'Position');
clrbarPos = get(clrbar,'Position');
set(clrbar,'Position',[pos1(1) .35*pos12(2) (pos12(1) + pos12(3) - pos1(1)) clrbarPos(4)])
set(get(clrbar,'xlabel'),'String', 'log_{10}[Observation Time (s)]','FontSize',12);
set(clrbar,'XAxisLocation','bottom')
set(get(clrbar,'XLabel'),'Position',get(get(clrbar,'XLabel'),'Position')+[0 3.5 0])
set(55,'Color',[1 1 1])

if saveplots
    if hex > 0
        export_fig HexagonalHERA_ObservationTimes.png  -nocrop -painters -r200
    else
        export_fig SquareHERA_ObservationTimes.png  -nocrop -painters -r200
    end
end

%%

figure(56); close(56); hfig = figure(56); clf;
set(56,'Position',[2059          89         678         859]);
handle = tight_subplot(4,3,[.07 .07],[.15 .05],[.06 .01]);
for f = 1:12
    axes(handle(f));
    PSF = fftshift(fft2(fftshift(observationTimesList{f})));
    PSF = PSF ./ max(max(PSF));
    imagesc(thetaLabelList{f},thetaLabelList{f},log10(abs(PSF)));
    colormap jet
    xlabel('\theta_x (degrees)'); ylabel('\theta_y (degrees)');
    title([num2str((80+10*f)) ' MHz'],'FontSize',14)
    caxis([-3 0])
end

axes(handle(1));
pos1 = get(handle(1),'Position');
clrbar = colorbar('location','southoutside');
set(handle(1),'Position',pos1);

pos12 = get(handle(12),'Position');
clrbarPos = get(clrbar,'Position');
set(clrbar,'Position',[pos1(1) .35*pos12(2) (pos12(1) + pos12(3) - pos1(1)) clrbarPos(4)])
set(get(clrbar,'xlabel'),'String', 'log_{10}[PSF]','FontSize',12);
set(clrbar,'XAxisLocation','bottom')
set(get(clrbar,'XLabel'),'Position',get(get(clrbar,'XLabel'),'Position')+[0 3.5 0])
set(56,'Color',[1 1 1])

if saveplots
    if hex > 0
        export_fig HexagonalHERA_PSF.png  -nocrop -painters -r200
    else
        export_fig SquareHERA_PSF.png  -nocrop -painters -r200
    end
end

