close all; clear all;

%% Settings

filename = 'PAPER112.dat';
xSpacing = 14.5;
ySpacing = 3.7;

%% Generate Array
array = [];
for x = 1:16
    for y = 1:7
        array = [array; (x-1)*xSpacing (y-1)*ySpacing];
    end
end

%% Print to File

fileID = fopen(filename,'w');
fprintf(fileID,'%u\n',length(array))
for (n = 1:length(array))
    fprintf(fileID,'%f %f\n',array(n,1),array(n,2));
end
fclose(fileID);