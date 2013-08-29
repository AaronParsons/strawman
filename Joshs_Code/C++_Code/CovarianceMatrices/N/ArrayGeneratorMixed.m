close all; clear all;

%% Settings

filename = 'Mixed_256.dat';
xSpacing = 4;
ySpacing = 4;

%% Generate Array
array = [];
for x = 1:31
    for y = 1:30
        if ((mod(x,3) == 1 && mod(y,3) == 1) || (x > 9 && x < 23 && y > 9 && y < 23))
            array = [array; (x-1)*xSpacing (y-1)*ySpacing];
        end
    end
end

array = [array; 120 120; 0 120];

plot(array(:,1),array(:,2),'.')
length(array)
%% Print to File

fileID = fopen(filename,'w');
fprintf(fileID,'%u\n',length(array))
for (n = 1:length(array))
    fprintf(fileID,'%f %f\n',array(n,1),array(n,2));
end
fclose(fileID);