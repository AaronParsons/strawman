%% POWER SPECTRUM SENSITIVITY ANALYSIS FOR HERA

clear all; close all;

ArrayFilenames = {'/Users/jsdillon/Desktop/HERA/HERA_Runs/HERA37',
    '/Users/jsdillon/Desktop/HERA/HERA_Runs/HERA127',
    '/Users/jsdillon/Desktop/HERA/HERA_Runs/HERA331',
    '/Users/jsdillon/Desktop/HERA/HERA_Runs/HERA547'};

nArrays = length(ArrayFilenames);

%% Settings
MproptoDelta = 0;
invVarWeightingSpherical = 1;
fields = 9;
exciseWedge = 1;
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

%% Loop over all arrays and foreground models;

vertErrorList = cell(3,nArrays);
kSphereBinCenterList = cell(3,nArrays);
leftHorizErrList = cell(3,nArrays);
rightHorizErrList = cell(3,nArrays);
WsphereList = cell(3,nArrays);

for a = 1:nArrays
    cd(ArrayFilenames{a});
    
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
    
    covP = M*fisher*M';
    
    %% Loop over Foreground models
    for foreModel = 1:3
        
        %% Set the foreground model
        if foreModel == 1
            wedgeEqualsHorizon = 1;
            buffer = .15;
            disp('Working on the conservative case...');
        elseif foreModel == 2
            wedgeEqualsHorizon = 1;
            buffer = .03;
            disp('Working on the fiducial case...');
        else
            wedgeEqualsHorizon = 0;
            buffer = .03;
            disp('Working on the optimistic case...');
        end
        
        %% Calculate Wedge
        PBangle = 13.15 * (f21cm/(zCenter + 1)/(150))^-1; %Degrees
        wedgeCoefficientPB = sin(PBangle*pi/180)*sqrt(OmegaM*(1+zCenter)^3+OmegaL)*DM/DH/(1+zCenter);
        wedgeCoefficientHoriz = sin(90*pi/180)*sqrt(OmegaM*(1+zCenter)^3+OmegaL)*DM/DH/(1+zCenter);
        if wedgeEqualsHorizon
            wedgeCoefficient = wedgeCoefficientHoriz;
        else
            wedgeCoefficient = wedgeCoefficientPB;
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
                    kValues(counter,1) = sqrt(kParaBinCenters(kPara)^2);% + kPerpBinCenters(kPerp)^2);
                    kValues(counter,2) = kPerp;
                    kValues(counter,3) = kPara;
                    counter = counter + 1;
                end
            end
        end
        kValues = sortrows(kValues);
        
        %Sets the number of kSphereBins
        binSize  = .063 / littleh;
        
        %Assigns bins in order
        assignment = zeros(kParaBins*kPerpBins,1);
        counter = 1;
        bin = 1;
        while true
            assignment(kValues(counter,3)+(kValues(counter,2)-1)*kParaBins) = kValues(counter,3);
            counter = counter + 1;
            if counter > size(kValues,1)
                break
            end
        end
        assignment = assignment - min(assignment(find(assignment))) + 1;
        
        %imagesc(reshape(assignment,kParaBins,kPerpBins))
        
        assignmentPlot = reshape(assignment,kParaBins,kPerpBins);
        
        assignmentMatrix = zeros(kParaBins+1,kPerpBins);
        
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
                    if (total >= percentile && left == 0)
                        leftHorizErr = [leftHorizErr; xi(i)];
                        left = 1;
                    end
                    if (total >= .5 && center == 0)
                        centerHorizErr = [centerHorizErr; xi(i)];
                        center = 1;
                    end
                    if (total >= (1-percentile) && right == 0)
                        rightHorizErr = [rightHorizErr; xi(i)];
                        right = 1;
                    end
                    if (~isnan(yi(i)))
                        total = total + yi(i)/norm;
                    end
                end
            end
            
        end
        kSphereBinCentersBackup = kSphereBinCenters;
        kSphereBinCenters = centerHorizErr;
        
        %% Lidz Theory
        
        lidzPS = importdata('/Users/jsdillon/Desktop/HERA/power_21cm_z7.32.dat');
        theoryCurve = [littleh*lidzPS(:,1)'; (lidzPS(:,1).^3/2/pi^2.*lidzPS(:,2)*(28 * sqrt( 9.5/10))^2)']';
        theoryInterpDeltaSq = interp1(theoryCurve(:,1),theoryCurve(:,2),kSphereBinCenters,'linear','extrap');
        
        
        %% Sample Variance
        
        independentBinCount = A'*reshape(cylindricalBinCount',kParaBins*kPerpBins,1)*fields;
        SampleVarianceDeltaSqError = (theoryInterpDeltaSq.*sqrt(2)./(independentBinCount).^.5);
        KPrefactor = kSphereBinCenters.^3 / (2*pi^2);
        vertError = 1000^2*((KPrefactor.*(diag(covPsphere)).^.5))/sqrt(fields);
        vertError = (vertError.^2 + SampleVarianceDeltaSqError.^2).^.5;
        
        vertErrorList{foreModel,a} = vertError;
        kSphereBinCenterList{foreModel,a} = kSphereBinCenters;
        leftHorizErrList{foreModel,a} = leftHorizErr;
        rightHorizErrList{foreModel,a} = rightHorizErr;
        WsphereListList{foreModel,a} = Wsphere;
        
    end
end

%% Plot Results: Foreground Comparison

arrayNames = {'HERA-37','HERA-127','HERA-331','HERA-547'};

figure(4); close(4); hfig = figure(4); clf;
set(4,'Position',[1816         852         743         703]);
ha = tight_subplot(2,2,[.0 .0],[.08 .08],[.08 .03]);


for a = 1:nArrays
    axes(ha(a));
    th = loglog(theoryCurve(:,1)/littleh, theoryCurve(:,2),'k','LineWidth',3);
    for m = 1:3
        kEdges = zeros(2*length(kSphereBinCenterList{m,a}),1);
        vertErrorLines = zeros(2*length(kSphereBinCenterList{m,a}),1);
        %kEdges(1:2:end) = leftHorizErrList{m,thisArray}/littleh;
        %kEdges(2:2:end) = rightHorizErrList{m,thisArray}/littleh;
        
        kEdges(1) = leftHorizErrList{m,a}(1)/littleh;
        for k = 1:length(kSphereBinCenterList{m,a})-1
            kEdges(2*k) = (kSphereBinCenterList{m,a}(k)+kSphereBinCenterList{m,a}(k+1))/2/littleh;
            kEdges(2*k+1) = (kSphereBinCenterList{m,a}(k)+kSphereBinCenterList{m,a}(k+1))/2/littleh;
        end
        kEdges(end) = kSphereBinCenterList{m,thisArray}(end)/littleh;
        kEdges = [kEdges(1); kEdges];
        
        
        vertErrorLines(1:2:end) = vertErrorList{m,a};
        vertErrorLines(2:2:end) = vertErrorList{m,a};
        vertErrorLines = [10^50; vertErrorLines];
        hold on;
        plotHandles{m} = loglog(kEdges, vertErrorLines,'b--');
        hold off
        %hold on; loglog(kSphereBinCenterList{m,thisArray}/littleh,vertErrorList{m,thisArray},'.'); hold off;
        set(gca,'XLim',[.03 2],'YLim',[.01 10^4],'XScale','log','YScale','log');
    end
    set(plotHandles{1},'Color','r','LineStyle','-.')
    set(plotHandles{2},'Color','k','LineStyle','-')
    if a == 1
        legend([plotHandles{1} plotHandles{2},plotHandles{3},th],'Conservative FGs', 'Fiducial FGs', 'Optimistic FGs','Lidz x_{HI} = 0.5','Location','SouthEast')
    end
    if a == 2 || a == 4
        set(gca,'YTickLabel',[]);
    else
        ylabel('\Delta^2(k)  (mK^2)')
    end
    if a == 1 || a == 2
        set(gca,'XTickLabel',[]);
    else
        xlabel('k (hMpc^{-1})')
    end    
        
end
changeFontSize(12)

for a = 1:nArrays
    axes(ha(a));
    text(.033,.02,arrayNames{a},'FontSize',20);
end

axes(ha(2)); 
hold on; txt = text(.03,10^4.5,['z = ' num2str(zRangeStart,3) ' to z = ' num2str(zRangeStop,3)],'HorizontalAlignment','center','FontSize',24);

set(gcf,'Color',[1 1 1])
export_fig /Users/jsdillon/Desktop/HERA/HERA_Foreground_Comp.pdf -r400 -nocrop -painters


%% Plot Results: Array Comparison

modelNames = {'Conservative FGs', 'Fiducial FGs', 'Optimistic FGs'};

figure(5); close(5); hfig = figure(5); clf;
set(5,'Position',[1816         -51         671        1029]);
ha = tight_subplot(3,1,[.0 .0],[.07 .05],[.1 .04]);


for m = 1:3
    axes(ha(m));
    th = loglog(theoryCurve(:,1)/littleh, theoryCurve(:,2),'k','LineWidth',3);
    for a = 1:nArrays
        
        kEdges = zeros(2*length(kSphereBinCenterList{m,a}),1);
        vertErrorLines = zeros(2*length(kSphereBinCenterList{m,a}),1);
        %kEdges(1:2:end) = leftHorizErrList{m,thisArray}/littleh;
        %kEdges(2:2:end) = rightHorizErrList{m,thisArray}/littleh;
        
        kEdges(1) = leftHorizErrList{m,a}(1)/littleh;
        for k = 1:length(kSphereBinCenterList{m,a})-1
            kEdges(2*k) = (kSphereBinCenterList{m,a}(k)+kSphereBinCenterList{m,a}(k+1))/2/littleh;
            kEdges(2*k+1) = (kSphereBinCenterList{m,a}(k)+kSphereBinCenterList{m,a}(k+1))/2/littleh;
        end
        kEdges(end) = kSphereBinCenterList{m,thisArray}(end)/littleh;
        kEdges = [kEdges(1); kEdges];
        
        
        vertErrorLines(1:2:end) = vertErrorList{m,a};
        vertErrorLines(2:2:end) = vertErrorList{m,a};
        vertErrorLines = [10^50; vertErrorLines];
        hold on;
        plotHandles{a} = loglog(kEdges, vertErrorLines,'k');
        hold off;
        %hold on; loglog(kSphereBinCenterList{m,thisArray}/littleh,vertErrorList{m,thisArray},'.'); hold off;
        set(gca,'XLim',[.03 2],'YLim',[.01 10^4],'XScale','log','YScale','log');
        
        
    end
    set(plotHandles{1},'Color','m','LineStyle',':')
    set(plotHandles{2},'Color','r','LineStyle','-.')
    set(plotHandles{3},'Color','b','LineStyle','--')
    if m == 3
        xlabel('k (hMpc^{-1})')
    else
        set(gca,'XTickLabel',[]);
    end
    ylabel('\Delta^2(k)  (mK^2)')
end
legend([plotHandles{1},plotHandles{2},plotHandles{3},plotHandles{4},th],'HERA-37','HERA-127','HERA-331','HERA-547','Lidz x_{HI} = 0.5','Location','SouthEast')


changeFontSize(16)

for m = 1:3
    axes(ha(m));
    text(.033,.02,modelNames{m},'FontSize',20);
end

axes(ha(1)); 
hold on; txt = text(sqrt(.03*2),10^4.5,['z = ' num2str(zRangeStart,3) ' to z = ' num2str(zRangeStop,3)],'HorizontalAlignment','center','FontSize',24);


set(gcf,'Color',[1 1 1])
export_fig /Users/jsdillon/Desktop/HERA/HERA_Array_Comp.pdf -r400 -nocrop -painters
