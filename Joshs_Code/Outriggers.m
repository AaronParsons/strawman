Rows = 24;
Cols = 24;
freq = 155*10^6;
Separation = 14;
Latitude = -30.7224 * 2*pi/360;
hourAngle = 0 * 2*pi/24;
timeSteps = 3;
c = 299792458;
totalTimeHours = 100;

saveplots = 0;
redundantOutriggers = 1;

hex = 50;
%hex = 0;
rotate = 0;

figure(1); clf

positions = [];



for row = hex-1 : -1 : -(hex-1)
    for col = 0:(2*hex-abs(row))-2
        xPos = ((-(2*hex-abs(row))+2)/2 + col)*sqrt(3)/2*Separation;
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
        xPos = ((-(2*hex-abs(row))+2)/2 + col)*sqrt(3)/2*Separation*outerPackingRatio;
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


figure(1);
plot(positions(:,1),positions(:,2),'ko','MarkerSize',4)
axis equal
xlabel('East-West Position (m)')
ylabel('North-South Position (m)')
plotSize = 1200;
set(gca,'XLim',[-plotSize plotSize],'YLim',[-plotSize plotSize])
set(1,'Color',[1 1 1])
changeFontSize(14);
if redundantOutriggers
    export_fig Outriggers.pdf  -nocrop
else
    export_fig Outriggers_Halved.pdf  -nocrop
end
    
%% Print to File

filename = 'Hex_HERA_With_Outriggers.txt';
fileID = fopen(filename,'w');
for (n = 1:length(positions))
    fprintf(fileID,'%f %f\n',positions(n,1),positions(n,2));
end
fclose(fileID);

%%

% 
% allBaselines = {};
% counter = 1;
% for i = 1:length(positions)
%     for j = 1:length(positions)
%         if i ~= j && ~(i > startOutriggers && j > startOutriggers)
%             allBaselines{counter,1} = (positions(i,1) - positions(j,1));
%             allBaselines{counter,2} =(positions(i,2) - positions(j,2));
%             counter = counter + 1;
%         end
%     end
% end
% uniqueBaselines = unique(cell2mat(allBaselines),'rows');
% 
% 

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
baselineKeys = keys(baselineMap);
baselineValues = values(baselineMap);
baselineList = [];
baselineSingleObs = [];
for k = 1:length(baselineKeys)
    commaPos = strfind(baselineKeys{k},',');
    deltaEast = str2num(baselineKeys{k}(1:(commaPos-1)));
    deltaNorth = str2num(baselineKeys{k}((commaPos+1):end));
    baselineList = [baselineList; deltaEast deltaNorth];
    if baselineValues{k} == 1 + redundantOutriggers
        baselineSingleObs = [baselineSingleObs; deltaEast deltaNorth];
    end
end

%%
figure(555); close(555); hfig = figure(555); clf;
set(555,'Position',[2007         210         867         723])
lambda = c/150000000;
plot(baselineList(:,1)/lambda, baselineList(:,2)/lambda,'o','MarkerSize',2,'MarkerFaceColor','b')
hold on; plot(baselineSingleObs(:,1)/lambda, baselineSingleObs(:,2)/lambda,'ro','MarkerSize',2,'MarkerFaceColor','r')
%plot(uniqueBaselines(:,1)/lambda, uniqueBaselines(:,2)/lambda,'.')
title(['Baseline Configuration at ' num2str((80+10*7)) ' MHz'],'FontSize',14)
xlabel('u'); ylabel('v');

for r = 100:100:600
    circle = [];
    for theta = 0:.01:2*pi
        circle = [circle; r*cos(theta) r*sin(theta)];
    end
    hold on; plot(circle(:,1),circle(:,2),'k--','LineWidth',2); hold off
    text(r,0,[num2str(r) '\lambda'],'BackgroundColor','w')
end

changeFontSize(14);
axis square
legend('Redundant Baselines','Non-Redundant Baselines')
set(555,'Color',[1 1 1])
export_fig Baselines.png  -nocrop -painters -r200

