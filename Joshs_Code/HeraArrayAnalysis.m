%% POWER SPECTRUM ANALYSIS FOR HERA

clear all; close all;
cd('/Users/jsdillon/Desktop/HERA/HERA_Runs/HERA547')

%% Settings
MproptoDelta = 0;
invVarWeightingSpherical = 1;
sampleVarianceOn = 1;

fields = 9;
exciseWedge = 1;
wedgeEqualsHorizon = 1;
buffer = .15;
buffer = .03;
nLowKPerpBinsToRemove = 0;
nHighKPerpBinsToRemove = 0;
nLowKParaBinsToRemove = 0;
folder = 'QuadraticEstimators/';
colorlimitsLog10K = [-6 1];
colorlimitsFisher = [-5 2];
colorlimitsPk = [-2 3];
SNRColors = [-2 2];
PowerSpectrumXLim = [10^-3.3 10^-.9];
PowerSpectrumYLim = [10^-2.5 10^.5];
SphereXLim = [10^-2 10^1];
SphereYLim = [10^-4 10^6];

%% Load in info about the data cube

disp(['Now working on ' pwd]);
fid = fopen([pwd '/QuadraticEstimators/cubeParameters.txt'],'r');
cubeParameters=textscan(fid,'%s',100,'delimiter','\n');


fStart = str2num(cubeParameters{1}{17});
fLength = str2num(cubeParameters{1}{14});
H0 = 70900;
littleh = H0/100000;
OmegaM = .27;
c = 299792000;
OmegaL = .73;
deltaRedshift = .00001;
f21cm = 1420.41;
comovingDist = 0;

zRight = f21cm/fStart - 1;
zLeft = zRight - deltaRedshift;
while (1) %Integrating backwards in z
    comovingDist = comovingDist + c/H0*((1.0/sqrt(OmegaM*(1+zLeft)^3 + OmegaL) + 4.0/sqrt(OmegaM*(1+(zLeft + zRight)/2)^3 + OmegaL) + 1.0/sqrt(OmegaM*(1+zRight)^3+OmegaL))*deltaRedshift/6);
    if (comovingDist >= fLength)
        break;
    end
    zRight = zLeft;
    zLeft = zLeft - deltaRedshift;
end
fL = f21cm/(zRight + 1);

zRangeStop = f21cm/fStart - 1;
zRangeStart = zRight;

DH = c/H0;
DM = 0;
zLeft = 0;
zRight = zLeft + deltaRedshift;
zCenter = (zRangeStop +zRangeStart)/2;
while (1)
    DM = DM + c/H0*((1.0/sqrt(OmegaM*(1+zLeft)^3 + OmegaL) + 4.0/sqrt(OmegaM*(1+(zLeft + zRight)/2)^3 + OmegaL) + 1.0/sqrt(OmegaM*(1+zRight)^3+OmegaL))*deltaRedshift/6);
    if (zRight > zCenter)
        break
    end
    zLeft = zRight;
    zRight = zLeft + deltaRedshift;
end


%% Load in Power Spectrum Estimates and compute Fisher Matrix

load QuadraticEstimators/fisher.mat;
load QuadraticEstimators/kParaBinCenters.dat;
load QuadraticEstimators/kPerpBinCenters.dat;
cylindricalBinCount = load('QuadraticEstimators/BinSampleCounts.dat');
%fisher = fisher*4;

%disp('Multiplying Fisher matrix by 4...')

%kParaBinCenters = kParaBinCenters / littleh;
%kPerpBinCenters = kPeBinCenters / littleh;


kPerpBins = length(kPerpBinCenters);
kParaBins = length(kParaBinCenters);

for n=1:nHighKPerpBinsToRemove
    fisher = fisher(1:kParaBins*kPerpBins - kParaBins, 1:kParaBins*kPerpBins - kParaBins);
    kPerpBinCenters = kPerpBinCenters(1:length(kPerpBinCenters)-1);
    bias = bias(1:(kPerpBins-1)*kParaBins);
    cylindricalBinCount = cylindricalBinCount(1:kPerpBins-1,:);
    kPerpBins = kPerpBins - 1;
end

while (max(max(abs(fisher))) / min(min(abs(fisher))) > 1e12)
    fisher = fisher(1:kParaBins*kPerpBins - kParaBins, 1:kParaBins*kPerpBins - kParaBins);
    kPerpBinCenters = kPerpBinCenters(1:length(kPerpBinCenters)-1);
    bias = bias(1:(kPerpBins-1)*kParaBins);
    cylindricalBinCount = cylindricalBinCount(1:kPerpBins-1,:);
    kPerpBins = kPerpBins - 1;
end

for n=1:nLowKPerpBinsToRemove
    fisher = fisher(kParaBins + 1:end,kParaBins + 1:end);
    kPerpBinCenters(1) = [];
    bias(1:kParaBins) = [];
    cylindricalBinCount(1,:) = [];
    kPerpBins = kPerpBins - 1;
end

for n=1:nLowKParaBinsToRemove
    kParaBinCenters(1) = [];
    fisher(1:kParaBins:end,:)=[];
    fisher(:,1:kParaBins:end)=[];
    bias(1:kParaBins:end)=[];
    cylindricalBinCount(end,:) = [];
    kParaBins = kParaBins -1;
end

%% Down-Sample Fisher Matrix

downPara = 1;
downPerp = 1;

newKParaBins = floor(kParaBins/downPara);
newKParaBinCenters = zeros(newKParaBins,1);
for l=1:newKParaBins
    newKParaBinCenters(l) = mean(kParaBinCenters((l-1)*downPara+1:l*downPara));
end
newKPerpBins = floor(kPerpBins/downPerp);
newKPerpBinCenters = zeros(newKPerpBins,1);
for m=1:newKPerpBins
    newKPerpBinCenters(m) = mean(kPerpBinCenters((m-1)*downPerp+1:m*downPerp));
end

newFisher = zeros(newKParaBins*newKPerpBins);
for mRow = 1:newKPerpBins*downPerp
    for lRow = 1:newKParaBins*downPara
        for mCol = 1:newKPerpBins*downPerp
            for lCol = 1:newKParaBins*downPara
                nRow = (mRow-1)*kParaBins + lRow;
                nCol = (mCol-1)*kParaBins + lCol;
                nDownRow = (ceil(mRow/downPerp)-1)*newKParaBins + ceil(lRow/downPara);
                nDownCol = (ceil(mCol/downPerp)-1)*newKParaBins + ceil(lCol/downPara);
                newFisher(nDownCol,nDownRow) = newFisher(nDownCol,nDownRow) + fisher(nCol,nRow);
            end
        end
    end
end

kParaBins = newKParaBins;
kParaBinCenters = newKParaBinCenters;
kPerpBins = newKPerpBins;
kPerpBinCenters = newKPerpBinCenters;
fisher = newFisher;


%% Set up results folder

currentDirectory = pwd;
[upperPath, ArrayName, ~] = fileparts(currentDirectory);

if MproptoDelta
    resultsFolderName = ['Corr Results ' ArrayName];
else
    resultsFolderName = ['Decorr Results ' ArrayName];
end

if (exist(resultsFolderName) == 7)
    rmdir(resultsFolderName,'s');
end
mkdir(resultsFolderName);


%% Fisher Diagonal

yaxislabels = [kParaBinCenters; kParaBinCenters(kParaBins)*2-kParaBinCenters(kParaBins-1)];
xaxislabels = kPerpBinCenters;
for n = length(xaxislabels):-1:2
    xaxislabels(n) = (xaxislabels(n) + xaxislabels(n-1))/2;
end
for n = length(yaxislabels):-1:2
    yaxislabels(n) = (yaxislabels(n) + yaxislabels(n-1))/2;
end
xaxislabels(1) = 2*xaxislabels(1)/3;
yaxislabels(1) = 2*yaxislabels(1)/3;


fisherDiag = zeros(kParaBins+1,kPerpBins);
fisherDiag(1:kParaBins, 1:kPerpBins) =  reshape(diag(fisher),kParaBins,kPerpBins);
fisherDiag(:,kPerpBins) = 0;
figure(3); close(3); hfig = figure(3); clf;
colormap jet;
psuedocolorpolot = pcolor(xaxislabels,yaxislabels,log10(fisherDiag)); %shading flat;
set(psuedocolorpolot,'LineWidth',.1);
clrbar = colorbar('peer',gca);
set(get(clrbar,'ylabel'),'String', 'log_{10}[diag(F)]  (Mpc^{-6} K^{-4})','FontSize',14);
locate = get(clrbar,'ylabel');
pos = get(locate,'Position');
set(locate,'Position',[pos(1)+1.5 pos(2) pos(3)])
caxis(colorlimitsFisher);
set(gca,'YScale','log','XScale','log');
set(gca,'XLim',PowerSpectrumXLim,'YLim',PowerSpectrumYLim);
xlabel('k_{\perp}  (Mpc^{-1})','FontSize',14);
ylabel('k_{||}  (Mpc^{-1})','FontSize',14);
title({[ArrayName ': Fisher Diagonal'],['z = ' num2str(zRangeStart,3) ' to z = ' num2str(zRangeStop,3)]});

changeFontSize(12);
set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
print(gcf, ['-r' num2str(400)], [resultsFolderName '/FisherDiagPlot'], ['-d' 'png']);


%% Calcualte M and W

sqrtF = sqrtm(fisher);
W = sqrtF;
invSqrtF = inv(sqrtm(fisher));
M = zeros(size(fisher));
sumSqrtF = zeros(kParaBins*kPerpBins,1);
sumF = zeros(kParaBins*kPerpBins,1);
for mCol = 1:kPerpBins
    for lCol = 1:kParaBins
        nCol = (mCol-1)*kParaBins + lCol;
        sumSqrtF(nCol) = real(sum(sqrtF(nCol,:)));
        sumF(nCol) = sum(fisher(nCol,:));
    end
end

for mRow = 1:kPerpBins
    for lRow = 1:kParaBins
        for mCol = 1:kPerpBins
            for lCol = 1:kParaBins
                nRow = (mRow-1)*kParaBins + lRow;
                nCol = (mCol-1)*kParaBins + lCol;
                if MproptoDelta
                    if (nCol == nRow)
                        M(nCol,nRow) = 1/sumF(nCol);
                    end
                    W(nCol,nRow) = fisher(nCol,nRow)/sumF(nCol);
                else
                    M(nCol,nRow) = invSqrtF(nCol,nRow)/sumSqrtF(nCol);
                    W(nCol,nRow) = sqrtF(nCol,nRow)/sumSqrtF(nCol);
                end
            end
        end
    end
end


%% Vertical Error Bars
covP = M*fisher*M';

% deltaPMatrix = zeros(kParaBins+1,kPerpBins);
% figure(5); close(5); hfig = figure(5); clf;
% colormap jet;
% n = 1;
% for m = 1:kPerpBins
%     for l = 1:kParaBins
%         deltaPMatrix(l,m) = abs((covP(n,n)))^.5;
%         n = n+1;
%     end
% end
% 
% psuedocolorpolot = pcolor(xaxislabels,yaxislabels,log10(deltaPMatrix)); %shading flat;
% set(psuedocolorpolot,'LineWidth',.1);
% caxis(colorlimitsPk);
% clrbar = colorbar('peer',gca);
% set(get(clrbar,'ylabel'),'String', 'log_{10}[Error on P(k)]  (K^2 cMpc^3)','FontSize',14);
% locate = get(clrbar,'ylabel');
% pos = get(locate,'Position');
% set(locate,'Position',[pos(1)+1.5 pos(2) pos(3)])
% 
% set(gca,'YScale','log','XScale','log');
% set(gca,'XLim',PowerSpectrumXLim,'YLim',PowerSpectrumYLim);
% xlabel('k_{\perp}  (Mpc^{-1})','FontSize',14);
% ylabel('k_{||}  (Mpc^{-1})','FontSize',14);
% title({[ArrayName ': Vertical Error Bars on P(k)'],['z = ' num2str(zRangeStart,3) ' to z = ' num2str(zRangeStop,3)]});
% 
% changeFontSize(12);
% set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
% if MproptoDelta
%     print(gcf, ['-r' num2str(400)], [resultsFolderName '/Corr_PkError'], ['-d' 'png']);
% else
%     print(gcf, ['-r' num2str(400)], [resultsFolderName '/Decorr_PkError'], ['-d' 'png']);
% end





deltaTMatrix = zeros(kParaBins+1,kPerpBins);
figure(1005); close(1005); hfig = figure(1005); clf;
%set(1005,'Position',[1864          92        1053         838])
set(1005,'Position',[560   403   776   545])
colormap jet;
n = 1;
for m = 1:kPerpBins
    for l = 1:kParaBins
        deltaTMatrix(l,m) = 1000^2*(kParaBinCenters(l)^2 + kPerpBinCenters(m)^2)^(3/2)/(2*pi^2)*abs((covP(n,n)))^.5;
        n = n+1;
    end
end

psuedocolorpolot = pcolor(xaxislabels/littleh,yaxislabels/littleh,log10(deltaTMatrix)); %shading flat;
set(psuedocolorpolot,'LineWidth',.1);
caxis([-3 4]);
clrbar = colorbar('peer',gca);
set(get(clrbar,'ylabel'),'String', 'log_{10}[Error on \Delta^2(k)]  (mK^2)','FontSize',14);
locate = get(clrbar,'ylabel');
pos = get(locate,'Position');
set(locate,'Position',[pos(1)+1.5 pos(2) pos(3)])

set(gca,'YScale','log','XScale','log');
set(gca,'XLim',[10^-3 10^-.7],'YLim',[.009 3]);
xlabel('k_{\perp}  (hMpc^{-1})','FontSize',14);
ylabel('k_{||}  (hMpc^{-1})','FontSize',14);
title({[ArrayName ' Simulation: Noise Vertical Error Bars on \Delta^2(k)'],['z = ' num2str(zRangeStart,3) ' to z = ' num2str(zRangeStop,3)]});


PBangle = 13.15 * (f21cm/(zCenter + 1)/(150))^-1; %Degrees
wedgeCoefficientPB = sin(PBangle*pi/180)*sqrt(OmegaM*(1+zCenter)^3+OmegaL)*DM/DH/(1+zCenter);
wedgeCoefficientHoriz = sin(90*pi/180)*sqrt(OmegaM*(1+zCenter)^3+OmegaL)*DM/DH/(1+zCenter);
if wedgeEqualsHorizon
    wedgeCoefficient = wedgeCoefficientHoriz;
else
    wedgeCoefficient = wedgeCoefficientPB;
end



hold on; plot([.001:.001:1], wedgeCoefficientHoriz*([.001:.001:1]) + .03,'w','LineWidth',4); hold off
text(.002,.07,'Horizon Wedge with .03 hMpc^{-1} Buffer','BackgroundColor',[1 1 1])
hold on; plot([.001:.001:1], wedgeCoefficientHoriz*([.001:.001:1]) + .03,'k--','LineWidth',2); hold off


hold on; plot([.001:.001:1], wedgeCoefficientPB*([.001:.001:1]) + .03,'w','LineWidth',4); hold off
text(.01,.02,'Primary Beam Null Wedge with .03 hMpc^{-1} Buffer','BackgroundColor',[1 1 1])
hold on; plot([.001:.001:1], wedgeCoefficientPB*([.001:.001:1]) + .03,'k--','LineWidth',2); hold off


hold on; plot([.001:.001:1], wedgeCoefficientHoriz*([.001:.001:1]) + .15,'w','LineWidth',4); hold off
text(.0015,.21,'Horizon Wedge with .15 hMpc^{-1} Buffer','BackgroundColor',[1 1 1])
hold on; plot([.001:.001:1], wedgeCoefficientHoriz*([.001:.001:1]) + .15,'k--','LineWidth',2); hold off


changeFontSize(14);
set(1005,'Color',[1 1 1])
thisFolder = pwd;
cd(resultsFolderName);
if MproptoDelta
    export_fig Corr_DeltaError.png -nocrop -painters -r200
else
    export_fig CylDeltaSqError.png -nocrop -painters -r200;
end
cd(thisFolder);



%% Show example window functions
if 0
    list = [];
    i = 2;
    j = 2;
    while i < kParaBins -1
        while j < kPerpBins -3
            list = [list; i j];
            j = round(1.5*j);
        end
        j = 2;
        i = round(1.5*i);
    end

    figure(2); clf;
    set(2,'Position',[762 350 660 638])

    interpFactor = 7;
    interpolatedKperp = interp(kPerpBinCenters,interpFactor,2,.5);
    interpolatedKpara = interp(kParaBinCenters,interpFactor,2,.5);

    for (i = 1:length(list))
        kPara = list(i,1);
        kPerp = list(i,2);
        interpLevel = 3;
        window = W(kPara + (kPerp-1)*kParaBins,:);
        window = reshape(window,kParaBins,kPerpBins);
        mask = (window > .05*max(max(window)));
        window = window.*mask;
        hold on;
        F = [.025 .05 .025; .05 .8 .05; .025 .05 .025];
        windowPlot = interp2(kPerpBinCenters,kParaBinCenters,window,interpolatedKperp,interpolatedKpara');
        contour(interpolatedKperp,interpolatedKpara,windowPlot,5); hold off;
        set(gca,'YScale','log','XScale','log');
        colormap jet;
        caxis([0 1]);
        windowColorbar = colorbar;
    end

    set(get(windowColorbar,'ylabel'),'String', 'Unitless Ratio');
    xlabel('k_{\perp}  (Mpc^{-1})','FontSize',14);
    ylabel('k_{||}  (Mpc^{-1})','FontSize',14);
    set(gca,'XLim',PowerSpectrumXLim,'YLim',PowerSpectrumYLim);
    grid minor;
    title({[ArrayName ': Example Window Functions'], ['z = ' num2str(zRangeStart,3) ' to z = ' num2str(zRangeStop,3)]});
    changeFontSize(12);
    set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
    if MproptoDelta
        saveas(gcf, [resultsFolderName '/Corr_Windows'], 'pdf') %Save figure
    else
        saveas(gcf, [resultsFolderName '/Decorr_Windows'], 'pdf') %Save figure
    end
end

%% Generate Rebinning Matrix

kSphereBins = kParaBins;

%Get kValues into a sorted list (excluding wedge, if necessary)
kParaBinCentersExtended = [kParaBinCenters; 2*kParaBinCenters(kParaBins) - kParaBinCenters(kParaBins-1)];
kValues = [];
counter = 1;
kValuesVector = [];
for kPara = 1:kParaBins
    for kPerp = 1:kPerpBins
        kValuesVector(kPara + (kPerp-1)*kParaBins) = sqrt(kParaBinCenters(kPara)^2 + kPerpBinCenters(kPerp)^2);
        if ((kParaBinCenters(kPara) > (wedgeCoefficient * (kPerpBinCenters(kPerp)) + buffer * littleh ) || ~exciseWedge) && kPerpBinCenters(kPerp) < 10^10 && sqrt(kParaBinCenters(kPara)^2 + kPerpBinCenters(kPerp)^2) > 0)
            kValues(counter,1) = sqrt(kParaBinCenters(kPara)^2 + kPerpBinCenters(kPerp)^2);
            kValues(counter,2) = kPerp;
            kValues(counter,3) = kPara;
            counter = counter + 1;
        end
    end
end
kValues = sortrows(kValues);

deltaKPara = (kParaBinCenters(2) - kParaBinCenters(1))*2;

%Assigns bins in order
assignment = zeros(kParaBins*kPerpBins,1);
for counter = 1:length(kValues)
    for bin = 2:kParaBins
        if kValues(counter,1) < kParaBinCenters(bin)
            assignment(kValues(counter,3)+(kValues(counter,2)-1)*kParaBins) = bin-1;
            break
        elseif bin == kParaBins
            assignment(kValues(counter,3)+(kValues(counter,2)-1)*kParaBins) = kParaBins;
        end
    end
end

assignment = assignment - (min((assignment == 0) * max(assignment) + assignment) - 1);
if exciseWedge
    assignment = assignment .* (assignment ~= min(assignment));
end

%Create rebinning matrix
kSphereBins = bin;
rebinningMatrix = zeros(kSphereBins,kParaBins*kPerpBins);
for kParaBin = 1:kParaBins
    for kPerpBin = 1:kPerpBins
        for sphereBin = 1:kSphereBins
            cylBin = kParaBin+(kPerpBin-1)*kParaBins;
            if assignment(kParaBin + (kPerpBin-1)*kParaBins) == sphereBin
                rebinningMatrix(sphereBin,cylBin) = 1;
            end
        end
    end
end
A = rebinningMatrix';



% kSphereBins = kParaBins;
% 
% 
% %Get kValues into a sorted list (excluding wedge, if necessary)
% kParaBinCentersExtended = [kParaBinCenters; 2*kParaBinCenters(kParaBins) - kParaBinCenters(kParaBins-1)];
% kValues = [];
% counter = 1;
% kValuesVector = [];
% for kPara = 1:kParaBins
%     for kPerp = 1:kPerpBins
%         kValuesVector(kPara + (kPerp-1)*kParaBins) = sqrt(kParaBinCenters(kPara)^2 + kPerpBinCenters(kPerp)^2);
%         if ((kParaBinCenters(kPara) > (wedgeCoefficient * (kPerpBinCenters(kPerp)) + buffer * littleh ) || ~exciseWedge) && kPerpBinCenters(kPerp) < 10^10 && sqrt(kParaBinCenters(kPara)^2 + kPerpBinCenters(kPerp)^2) > 0)
%             kValues(counter,1) = sqrt(kParaBinCenters(kPara)^2);% + kPerpBinCenters(kPerp)^2);
%             kValues(counter,2) = kPerp;
%             kValues(counter,3) = kPara;
%             counter = counter + 1;
%         end
%     end
% end
% kValues = sortrows(kValues);
% 
% %Sets the number of kSphereBins
% binSize  = .063 / littleh;
% 
% %Assigns bins in order
% assignment = zeros(kParaBins*kPerpBins,1);
% counter = 1;
% bin = 1;
% while true
%         assignment(kValues(counter,3)+(kValues(counter,2)-1)*kParaBins) = kValues(counter,3);
%         counter = counter + 1;
%     if counter > size(kValues,1)
%         break
%     end
% end
% assignment = assignment - min(assignment(find(assignment))) + 1;
% 
% %imagesc(reshape(assignment,kParaBins,kPerpBins))
% 
assignmentPlot = reshape(assignment,kParaBins,kPerpBins);

assignmentMatrix = zeros(kParaBins+1,kPerpBins);
figure(105); close(105); hfig = figure(105); clf;
%set(105,'Position',[1864          92        1053         838])
set(105,'Position',[560   403   776   545])
customLines = lines;
customLines(1,:) = [0 0 0];
colormap(customLines)
n = 1;
for m = 1:kPerpBins
    for l = 1:kParaBins
        assignmentMatrix(l,m) = assignmentPlot(l,m);
        n = n+1;
    end
end

psuedocolorpolot = pcolor(xaxislabels/littleh,yaxislabels/littleh,assignmentMatrix); %shading flat;
set(psuedocolorpolot,'LineWidth',.1);
caxis([0 50]);
%clrbar = colorbar('peer',gca);
%set(get(clrbar,'ylabel'),'String', 'Bin Number','FontSize',14);
%locate = get(clrbar,'ylabel');
%pos = get(locate,'Position');
%set(locate,'Position',[pos(1)+1.5 pos(2) pos(3)])

set(gca,'YScale','log','XScale','log');
set(gca,'XLim',[10^-3 10^-.7],'YLim',[.009 3]);
xlabel('k_{\perp}  (hMpc^{-1})','FontSize',14);
ylabel('k_{||}  (hMpc^{-1})','FontSize',14);
title({[ArrayName ' Simulation: Binning Diagram'],['z = ' num2str(zRangeStart,3) ' to z = ' num2str(zRangeStop,3)]});


PBangle = 9.8 * (f21cm/(zCenter + 1)/(150))^-1; %Degrees
wedgeCoefficientPB = sin(PBangle*pi/180)*sqrt(OmegaM*(1+zCenter)^3+OmegaL)*DM/DH/(1+zCenter);
wedgeCoefficientHoriz = sin(90*pi/180)*sqrt(OmegaM*(1+zCenter)^3+OmegaL)*DM/DH/(1+zCenter);
if wedgeEqualsHorizon
    wedgeCoefficient = wedgeCoefficientHoriz;
else
    wedgeCoefficient = wedgeCoefficientPB;
end



hold on; plot([.001:.001:1], wedgeCoefficientHoriz*([.001:.001:1]) + .03,'w','LineWidth',4); hold off
text(.002,.07,'Horizon Wedge with .03 hMpc^{-1} Buffer','BackgroundColor',[1 1 1])
hold on; plot([.001:.001:1], wedgeCoefficientHoriz*([.001:.001:1]) + .03,'k--','LineWidth',2); hold off


hold on; plot([.001:.001:1], wedgeCoefficientPB*([.001:.001:1]) + .03,'w','LineWidth',4); hold off
text(.01,.02,'Primary Beam Null Wedge with .03 hMpc^{-1} Buffer','BackgroundColor',[1 1 1])
hold on; plot([.001:.001:1], wedgeCoefficientPB*([.001:.001:1]) + .03,'k--','LineWidth',2); hold off


hold on; plot([.001:.001:1], wedgeCoefficientHoriz*([.001:.001:1]) + .15,'w','LineWidth',4); hold off
text(.0015,.21,'Horizon Wedge with .15 hMpc^{-1} Buffer','BackgroundColor',[1 1 1])
hold on; plot([.001:.001:1], wedgeCoefficientHoriz*([.001:.001:1]) + .15,'k--','LineWidth',2); hold off



changeFontSize(14);
set(105,'Color',[1 1 1])
thisFolder = pwd;
cd(resultsFolderName);
if MproptoDelta
    export_fig Binning.png -nocrop -painters -r200
else
    export_fig Binning.png -nocrop -painters -r200;
end
cd(thisFolder);


%%

% %Sets the number of kSphereBins
% binExpansionFactor = 1.8;
% 
% %Assigns bins in order
% assignment = zeros(kParaBins*kPerpBins,1);
% counter = 1;
% bin = 1;
% while true
%     for i = 1:(6+ceil(binExpansionFactor^bin))
%         assignment(kValues(counter,3)+(kValues(counter,2)-1)*kParaBins) = bin;
%         counter = counter + 1;
%         if counter > size(kValues,1)
%             break
%         end
%     end
%     if counter > size(kValues,1)
%         break
%     end
%     bin = bin + 1;
% end
% imagesc(reshape(assignment,kParaBins,kPerpBins))




%Create rebinning matrix
kSphereBins = max(assignment);
rebinningMatrix = zeros(kSphereBins,kParaBins*kPerpBins);
for kParaBin = 1:kParaBins
    for kPerpBin = 1:kPerpBins
        for sphereBin = 1:kSphereBins
            cylBin = kParaBin+(kPerpBin-1)*kParaBins;
            if assignment(kParaBin + (kPerpBin-1)*kParaBins) == sphereBin
                rebinningMatrix(sphereBin,cylBin) = 1;
            end
        end
    end
end
A = rebinningMatrix';

%% Lidz Theory

lidzPS = importdata('/Users/jsdillon/Desktop/HERA/power_21cm_z7.32.dat');
theoryCurve = [littleh*lidzPS(:,1)'; (lidzPS(:,1).^3/2/pi^2.*lidzPS(:,2)*(28 * sqrt( 9.5/10))^2)']';


%% Sample Variance
allSphereKValues = zeros(kParaBins*kPerpBins,1);
for kParaBin = 1:kParaBins
    for kPerpBin = 1:kPerpBins
        allSphereKValues(kParaBin + (kPerpBin-1)*kParaBins) = sqrt(kParaBinCenters(kParaBin)^2 + kPerpBinCenters(kPerpBin)^2);
    end
end

theoryForAllSphereKValues = interp1(theoryCurve(:,1),theoryCurve(:,2),allSphereKValues,'linear','extrap');
binCountVector = reshape(cylindricalBinCount',kParaBins*kPerpBins,1)/2;
cylindricalSampleVarianceDeltaSqError = theoryForAllSphereKValues.*sqrt(2)./binCountVector.^.5; %in mK^2
cylindricalSampleVariancePError = cylindricalSampleVarianceDeltaSqError ./ allSphereKValues.^3 * 2 * pi^2 / (1000)^2;

if sampleVarianceOn
    %covP = covP + W*diag(cylindricalSampleVariancePError.^2)*W';
    covP = covP + diag(cylindricalSampleVariancePError.^2);
    disp(['Not handling window functions properly for sample variance...I think...'])
end
    
 
% independentBinCount = A'*reshape(cylindricalBinCount',kParaBins*kPerpBins,1)*fields/2;
% SampleVarianceDeltaSqError = (theoryInterpDeltaSq.*sqrt(2)./(independentBinCount).^.5);
% 
% figure(321); clf;
% plot(kSphereBinCenters/littleh,independentBinCount,'.-')
% xlabel('k (hMpc^{-1})');
% ylabel('Number of independent bins');


%% Recalculate Band Powers, Window Functions, and Band Power Covariance

if invVarWeightingSpherical
    invVar = inv(covP);
    covPsphere = inv(A'*invVar*A);
    Wsphere = inv(A'*invVar*A)*A'*invVar*W*A;
    kSphereBinCenters = covPsphere*A'*invVar * kValuesVector';
else
    pHatSphereList = [];
    covPsphere = A'*covP*A;
    Wsphere = A'*W*A;
end

%% Adrian's N_indep Calculation

if invVarWeightingSpherical
    Nindep = (sum((covPsphere * A' * invVar).^2,2)).^-1;
    figure(1234); clf
    plot(kSphereBinCenters/littleh, Nindep,'.-')
    ylabel('N_{indep}'); xlabel('k (hMpc^{-1})');
end

% if invVarWeightingSpherical
%     b = (covPsphere * A' * invVar);
%     
%     figure(1235); clf
%     theorySq = (interp1(theoryCurve(:,1),theoryCurve(:,2),kSphereBinCenters,'linear','extrap') ./ kSphereBinCenters.^3 * 2 * pi^2 / (1000)^2).^2
%     plot(kSphereBinCenters/littleh,diag(b *  W*diag(cylindricalSampleVariancePError.^2)*W' * b') ./ theorySq)
%   
%     ylabel('N_{indep}'); xlabel('k (hMpc^{-1})');
% end

%% Calculate Horizontal Error Bars

leftHorizErr = [];
rightHorizErr = [];
centerHorizErr = [];

%percentile = 0.158655; %1 sigma
percentile = .2; % 20%-80%

stepSize = log10(kSphereBinCenters(end) - kSphereBinCenters(1))/10000;

for win = 1:kSphereBins
    xi = 10.^(log10(kSphereBinCenters(1)):stepSize:log10(kSphereBinCenters(end))) + 10e-15;
    yi = interp1(kSphereBinCenters,Wsphere(win,:),xi,'linear');
    ni = length(yi);
    norm = sum(yi);
    total = 0;
    left = 0;
    center = 0;
    right = 0;
    if win == 1
        centerHorizErr = [centerHorizErr; xi(1)];
        for i=1:ni
            if (total >= (3/4-percentile) && right == 0)
                rightHorizErr = [rightHorizErr; xi(i)];
                right = 1;
            end
            if (~isnan(yi(i)))
                total = total + yi(i)/norm;
            end
        end
        leftHorizErr = [leftHorizErr; max(2*centerHorizErr(end) - rightHorizErr(end),10e-15)];
    elseif win == kSphereBins
        centerHorizErr = [centerHorizErr; xi(end)];
        for i=1:ni
            if (total >= (1/4+percentile) && left == 0)
                leftHorizErr = [leftHorizErr; xi(i)];
                left = 1;
            end
            if (~isnan(yi(i)))
                total = total + yi(i)/norm;
            end
        end
        rightHorizErr = [rightHorizErr; max(2*centerHorizErr(end) - leftHorizErr(end),10e-15)];
    else
        for i=1:ni

            if (total >= .5 && center == 0)
                centerHorizErr = [centerHorizErr; xi(i)];
                center = 1;
            end
            if (total >= (1-percentile) && right == 0)
                rightHorizErr = [rightHorizErr; xi(i)];
                right = 1;
            end
            if (total >= percentile && left == 0)
                leftHorizErr = [leftHorizErr; max(xi(i),leftHorizErr(end))];
                left = 1;
            end
            if (~isnan(yi(i)))
                total = total + yi(i)/norm;
            end
        end
    end
    
end
kSphereBinCentersBackup = kSphereBinCenters;
%kSphereBinCenters = centerHorizErr;
disp('No longer using center of window functions as bin center')



%% Detections and Upper Limits

figure(22); close(22); hfig = figure(22); clf;
set(22,'Position',[560   256   804   692])
KPrefactor = kSphereBinCenters.^3 / (2*pi^2);
vertError = 1000^2*((KPrefactor.*(diag(covPsphere)).^.5))/sqrt(fields);
theoryInterpDeltaSq = interp1(theoryCurve(:,1),theoryCurve(:,2),kSphereBinCenters,'linear','extrap');
plotBottom = min(min(theoryCurve(:,2)),min(vertError/3));

leftBar = leftHorizErr/littleh;
rightBar = rightHorizErr/littleh;
bottomBar = max(theoryInterpDeltaSq - vertError,1e-50);
topBar = vertError + theoryInterpDeltaSq;


xHandleSize = .2*log10(SphereXLim(2)/SphereXLim(1))/20;
yHandleSize = 1.0*0.0048*log10(max(topBar)*1000/plotBottom);
h = ploterr(kSphereBinCenters/littleh, theoryInterpDeltaSq, {leftBar, rightBar}, {bottomBar,topBar},'ko','logx','logy','abshhy',yHandleSize,'abshhx',xHandleSize);
set(h(1),'MarkerSize',4,'MarkerFaceColor','k')


hold on; h2 = loglog(theoryCurve(:,1)/littleh,theoryCurve(:,2),'b--'); hold off;
set(gca,'XLim',[.02 1],'YLim',[1 10^2]);
xlabel('k (h Mpc^{-1})');
ylabel('\Delta^2(k)  (mK^2)');
title({[ArrayName ' Simulation of Errors on \Delta^2(k)'], ['z = ' num2str(zRangeStart,3) ' to z = ' num2str(zRangeStop,3)]});
set(gca,'YGrid','on','YMinorGrid','off');

for n = 1:length(kSphereBinCenters)
    text(kSphereBinCenters(n)/littleh,bottomBar(n)/1.1, num2str(theoryInterpDeltaSq(n)/vertError(n),3),'HorizontalAlignment','center')
end

%legend([h(1) h(3) h2],'\Delta^2(k) 1\sigma Errors','20%-80% Window Functions', 'Theory at z=8)', 'Location','NorthWest')
changeFontSize(14);
set(gca,'YMinorTick','on');
set(gca,'YTickMode','manual');

set(22,'Color',[1 1 1])
if MproptoDelta
    filename = [resultsFolderName '/Corr_Deltak.png'];
else
    filename = [resultsFolderName '/HERA_DeltaSqPrediction.png'];
end
export_fig(filename, '-nocrop')

sqrt(sum((theoryInterpDeltaSq./vertError).^2))
%% 1D Window Functions

figure(444); close(444); hfig = figure(444); clf;
semilogx(kSphereBinCentersBackup,Wsphere','.-','MarkerSize',16,'LineWidth',1.5)
set(gca,'XLim',SphereXLim,'YLim',[min(min(min(Wsphere))*2,0) 1]);
title({[ArrayName ': Window Functions on Spherical P(k)'], ['z = ' num2str(zRangeStart,3) ' to z = ' num2str(zRangeStop,3)]});
set(gca,'YGrid','on','YMinorGrid','off');
xlabel('k (Mpc^{-1})');
ylabel('Window Function (Unitless)');
changeFontSize(14);
set(gca,'YMinorTick','on');
set(gca,'YTickMode','manual');
set(gcf, 'PaperPosition', [0 0 6 6]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [6 6]); %Set the paper to have width 5 and height 5.
if MproptoDelta
    saveas(gcf, [resultsFolderName '/Corr_1D_Windows'], 'pdf') %Save figure
else
    saveas(gcf, [resultsFolderName '/Decorr_1D_Windows'], 'pdf') %Save figure
end

%%
wedgeEqualsHorizon
buffer
AllKBins = num2str(kSphereBinCenters/littleh)
DeltaSqError = num2str(vertError)
wedgeEqualsHorizon
buffer

%%

load('samp.txt')
load('samp-fg.txt')

figure(50)
plot(samp_fg(:,1),abs(samp_fg(:,2)),'.')


joshsCoords = [];
for kParaBin = 1:kParaBins
    for kPerpBin = 1:kPerpBins
        joshsCoords = [joshsCoords; (assignment(kParaBin + (kPerpBin-1)*kParaBins) ~= 0) * kPerpBinCenters(kPerpBin)/littleh (assignment(kParaBin + (kPerpBin-1)*kParaBins) ~= 0)*kParaBinCenters(kParaBin)/littleh];
    end
end
%joshsCoords(:,1) = joshsCoords(:,1).*(assignment ~= 0);
%joshsCoords(:,2) = joshsCoords(:,2).*(assignment ~= 0);

plot(joshsCoords(:,1),abs(joshsCoords(:,2)),'.')

hold on; plot(samp_fg(:,1),abs(samp_fg(:,2)),'r.'); hold off

hold on; plot([0 max(samp_fg(:,1))],[0 max(samp_fg(:,1))]*wedgeCoefficient+.15); hold off

xlabel('k_\perp (hMpc^{-1})'); ylabel('k_{||} (hMpc^{-1})')
changeFontSize(20)

saveas(gcf, [resultsFolderName '/JoshVsJonnieBinning'], 'pdf') %Save figure

legend('Josh','Jonnie')
