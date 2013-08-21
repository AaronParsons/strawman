close all; clear all;

%% Settings


Rows = 24;
Cols = 24;
filename = 'HERA37_HexPack_HexPerim_14m_Layout.dat';
hex = 4;
Separation = 14;


% Generate Array

figure(1); clf

array = [];
for row = hex-1 : -1 : -(hex-1)
    for col = 0:(2*hex-abs(row))-2
        xPos = ((-(2*hex-abs(row))+2)/2 + col)*Separation;
        yPos = row*Separation*sqrt(3)/2;
        radius = sqrt(xPos^2 + yPos^2);
        if radius < 175
            array = [array; xPos yPos];
        end
    end
end

length(array)

%array = array(1:Rows*Cols+3,:);

plot(array(:,1),array(:,2),'ko','MarkerSize',8)
axis equal
xlabel('East-West Position (m)')
ylabel('North-South Position (m)')
plotSize = 400;
set(gca,'XLim',[-plotSize plotSize],'YLim',[-plotSize plotSize])
set(1,'Color',[1 1 1])



%% Print to File

fileID = fopen(filename,'w');
fprintf(fileID,'%u\n',length(array));
for (n = 1:length(array))
    fprintf(fileID,'%f %f\n',array(n,1),array(n,2));
end
fclose(fileID);