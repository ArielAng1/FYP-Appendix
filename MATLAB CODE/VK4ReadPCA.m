%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   File displayer for Keyence laser conformal microscope (LCM)
%   File extension: vk4
%   
%	Adapted from: 
%   https://permalink.lanl.gov/object/tr?what=info:lanl-repo/lareport/LA-UR-18-20810
%   WTP 15/6/23
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WYKO: Read ASCII triplet file
clear all;
myfile = uigetfile('*.*');

% Define the filename and the number of header lines

% Open the file
fid = fopen(myfile, 'r');

% Skip the first four lines
for i = 1:4
    line = fgetl(fid);

end
% Initialize an empty structure to store the values
v = struct();

% Read the next 132 lines
while true
    % Get the next line from the file
    line = fgetl(fid);
    
    % Check if end of file is reached before 132 lines
    if line == -1
        break;
    end

    % Check if the line contains the text "Mult"
    if contains(line, 'Mult')
        %Break after reading the last line to get the multiplier
        % Split the line by tab delimiter
        splitLine = strsplit(line, '\t');
    
        % Check if the line has at least four columns
        if length(splitLine) >= 4
            varName = splitLine{1};  % First column: Variable name
            % Remove spaces from the variable name
            % Remove spaces from the variable name
            varName = regexprep(varName, '[^a-zA-Z0-9]', '');
            varValue = splitLine{4};  % Fourth column: Variable value initially as text)
            
            % Convert to number if possible
            numValue = str2double(varValue);
            
            if isnan(numValue)
                % If not a number, store as text
                v.(varName) = varValue;
            else
                % If it's a number, store as numeric
                v.(varName) = numValue;
            end
        end
        break; % Stop reading further lines
    end
    
    % Split the line by tab delimiter
    splitLine = strsplit(line, '\t');
    
    % Check if the line has at least four columns
    if length(splitLine) >= 4
        varName = splitLine{1};  % First column: Variable name
        % Remove spaces from the variable name
        % Remove spaces from the variable name
        varName = regexprep(varName, '[^a-zA-Z0-9]', '');
        varValue = splitLine{4};  % Fourth column: Variable value initially as text)
        
        % Convert to number if possible
        numValue = str2double(varValue);
        
        if isnan(numValue)
            % If not a number, store as text
            v.(varName) = varValue;
        else
            % If it's a number, store as numeric
            v.(varName) = numValue;
        end
    end
end

% Set image size
x = strsplit(v.XYSize);
xSize = str2double(x{1});
ySize = str2double(x{3});

l_per_pix(1) = v.Pixelsize;
l_per_pix(2) = v.Pixelsize;
wavelength = v.Wavelength;%
%changing from integer to floating point
scale = v.Mult;


%construct pixel vectors in µm
xPix = (1:ySize)*l_per_pix(1)*1000;
yPix = (1:xSize)*l_per_pix(2)*1000;


% Read the data
% The format specifier '%f %f %f' reads three floating-point numbers
% The '*' before '%f' tells MATLAB to ignore any missing value
data = textscan(fid, '%f %f %f', 'Delimiter', '\t', 'TreatAsEmpty', ' ');

% Close the file
fclose(fid);

% Convert the cell array to a matrix
X = data{1};
Y = data{2};
H = data{3};

heightVector = sort(H);
[f, xi] = ksdensity(heightVector);
H = H - xi(f == max(f));


H = [H; NaN(xSize*ySize-length(H), 1)];

myheight =reshape(H, [ySize, xSize])*wavelength/1000 /scale;

figure(4)
clf
imshow(rescale(myheight'/max(max(myheight))));

colormap('turbo')
view(-90,90)

%% Interpolate NaN values

% Create a copy of the original matrix to use for checking adjacent values
originalHeight = myheight;

% Get the size of the matrix
[numRows, numCols] = size(myheight);

% Iterate through each element in the matrix
for row = 1:numRows
    for col = 1:numCols
        if isnan(myheight(row, col))
            % Check for non-NaN values in the adjacent positions in the original matrix
            hasAdjacentValue = false;
            
            % Check the left value
            if col > 1 && ~isnan(originalHeight(row, col - 1))
                hasAdjacentValue = true;
            end
            
            % Check the right value
            if col < numCols && ~isnan(originalHeight(row, col + 1))
                hasAdjacentValue = true;
            end
            
            % Check the upper value
            if row > 1 && ~isnan(originalHeight(row - 1, col))
                hasAdjacentValue = true;
            end
            
            % Check the lower value
            if row < numRows && ~isnan(originalHeight(row + 1, col))
                hasAdjacentValue = true;
            end
            
            % If there is at least one non-NaN adjacent value, interpolate
            if hasAdjacentValue
                % Gather the indices and values of surrounding non-NaN points
                x = [];
                v = [];
                
                % Left adjacent value
                if col > 1 && ~isnan(originalHeight(row, col - 1))
                    x = [x, col - 1];
                    v = [v, originalHeight(row, col - 1)];
                end
                
                % Right adjacent value
                if col < numCols && ~isnan(originalHeight(row, col + 1))
                    x = [x, col + 1];
                    v = [v, originalHeight(row, col + 1)];
                end
                
                % Upper adjacent value
                if row > 1 && ~isnan(originalHeight(row - 1, col))
                    x = [x, numCols * (row - 2) + col];
                    v = [v, originalHeight(row - 1, col)];
                end
                
                % Lower adjacent value
                if row < numRows && ~isnan(originalHeight(row + 1, col))
                    x = [x, numCols * row + col];
                    v = [v, originalHeight(row + 1, col)];
                end
                
                % Interpolate based on adjacent values
                if ~isempty(x) && ~isempty(v)
                    % Handle interpolation for both row and column neighbors
                    if numel(x) == 1
                        % If only one valid neighbor, just copy its value
                        myheight(row, col) = v;
                    else
                        % Use linear interpolation for multiple neighbors
                        myheight(row, col) = mean(v);
                    end
                end
            end
        end
    end
end

% Erode to NaN
finalMatrix = myheight;

% Get the size of the matrix
[numRows, numCols] = size(finalMatrix);

% Create a logical matrix to mark elements that should be set to NaN
erodeMask = false(numRows, numCols);

% Iterate through each element in the matrix
for row = 1:numRows
    for col = 1:numCols
        if isnan(finalMatrix(row, col))
            % Mark adjacent elements for erosion
            if col > 1  % Left
                erodeMask(row, col - 1) = true;
            end
            if col < numCols  % Right
                erodeMask(row, col + 1) = true;
            end
            if row > 1  % Up
                erodeMask(row - 1, col) = true;
            end
            if row < numRows  % Down
                erodeMask(row + 1, col) = true;
            end
        end
    end
end

% Apply the erosion: set elements to NaN where the mask is true
myheight(erodeMask) = NaN;

figure(5)
clf
imshow(rescale(myheight'/max(max(myheight))));
colormap('jet')
view(-90,90)

%% PCA Leveling
[Xgrid, Ygrid] = meshgrid(1:size(myheight, 2), 1:size(myheight, 1));


valid = ~isnan(myheight);
X = Xgrid(valid);
Y = Ygrid(valid);
Z = myheight(valid);


points = [X, Y, Z];
meanVals = mean(points, 1);
centeredPoints = points - meanVals;


[~, ~, V] = svd(centeredPoints, 0);
normalVec = V(:, 3);  % The normal vector to the best-fit plane

% Reconstruct the best-fit plane
% z = -(nx*(x - x0) + ny*(y - y0))/nz + z0
nx = normalVec(1);
ny = normalVec(2);
nz = normalVec(3);

% Avoid division by zero
if abs(nz) < 1e-6
    error('Plane normal is nearly vertical. Cannot compute leveling plane.');
end

% Recreate the full plane using the normal vector
plane = -(nx * (Xgrid - meanVals(1)) + ny * (Ygrid - meanVals(2))) / nz + meanVals(3);

% Step 5: Subtract the plane to level the surface
leveledHeight = myheight - plane;

% Optional: visualize results
figure;
subplot(1,2,1);
imagesc(myheight);
title('Original Height Map');
axis image; colorbar;

subplot(1,2,2);
imagesc(leveledHeight);
title('Planar-Leveled Height Map (PCA)');
axis image; colorbar;

% Optional: output leveled data
myheight = leveledHeight;

%% Display 3D data

% % Vertical exaggeration
vertEx = 5; 
% colormap('turbo')
% figure(5)
% % Vectorize the height data
% heightVector = myheight(:);
% % Remove NaNs
% heightVector = heightVector(~isnan(heightVector));
% 
% % set limits to 1 and 99th percentile
% heightSort = sort(heightVector);
% heightLim = [heightSort(round(end/100)) heightSort(round(99*end/100))];
% 
% 
% %   Height image 
% h = surf(yPix, xPix, myheight,...
%     'EdgeColor','none');

% Downsample factor (adjust as needed)
ds = 1;
xPix_ds = xPix(1:ds:end, 1:ds:end);
yPix_ds = yPix(1:ds:end, 1:ds:end);
myheight_ds = myheight(1:ds:end, 1:ds:end);

% Use downsampled data in plotting
figure(5)
heightVector = myheight_ds(:);
heightVector = heightVector(~isnan(heightVector));
heightSort = sort(heightVector);
heightLim = [heightSort(round(end/100)) heightSort(round(99*end/100))];

h = surf(yPix_ds, xPix_ds, myheight_ds, 'EdgeColor', 'none');

camlight('right');
%camlight('left');
%camlight('headlight');
% camlight('headlight');
% camlight('headlight');
% camlight('headlight');
% camlight('headlight');
% camlight('headlight');


view([1 -1 1]);
    
set(gca, 'Projection','perspective')
axis tight
xlabel('µm')
ylabel('µm')
zlabel('µm')
daspect([1 1 1/vertEx]);
caxis(heightLim)
cb = colorbar
cb.Label.String = 'height/µm';


%% Display height data only
nHeight = 8192;
vertEx = 15;

colormap('hsv')
%construct pixel vectors in µm


figure(6)
clf
t = tiledlayout(3,6,'TileSpacing','Compact','Padding','Compact');
% Vectorize the height data
heightVector = myheight(:);
% Remove NaNs
heightVector = heightVector(~isnan(heightVector));
% set limits to 1 and 99th percentile
heightSort = sort(heightVector);
heightLim = [heightSort(round(end/100)) heightSort(round(99*end/100))];

xi = linspace(heightLim(1), heightLim(2), nHeight);

% Kernel density estimation
f = ksdensity(heightVector, xi);

myheight = myheight - xi(f == max(f));

% Find peaks in the kernel density estimate
[~, locs] = findpeaks(f, xi, 'MinPeakProminence', 0.01);
t4 =nexttile(13)
plot(xi, f, 'k', 'LineWidth', 1.5)
for i = 1:length(locs)
    xline(locs(i), 'k', num2str(locs(i), '%.2f'))
end 
xlabel('Height/µm')
ylabel ('Kernel pdf')
title('Height distribution of full dataset')
xlim(heightLim)

%   Height image 


t1 = nexttile(2,[2 5])
h = surf(yPix, xPix, myheight,...
    'EdgeColor','none');

% shading interp
% lightangle(-45,30)
 
% h.FaceLighting = 'gouraud';
% h.AmbientStrength = 0.5;
% h.DiffuseStrength = 0.8;
% h.SpecularStrength = 0.9;
% h.SpecularExponent = 25;
% h.BackFaceLighting = 'unlit';
camlight('right');
camlight('left');
camlight('headlight');
camlight('headlight');
camlight('headlight');
camlight('headlight');
camlight('headlight');
camlight('headlight');

view([0 0 1]);

set(gca, 'Projection','perspective')
axis tight
xlabel('µm')
ylabel('µm')
zlabel('µm')
daspect([1 1 1/vertEx]);
caxis(heightLim)
cb = colorbar
cb.Label.String = 'height/µm';
grid off

% direction
xX = 1;


% bands
profx = round([0.74*ySize, 0.75*ySize;...
               0.70*ySize, 0.71*ySize]);

profy = round([0.65*xSize, 0.66*xSize]);


cmapx = cool(size(profx,1))
for i = 1:size(profx,1)
yline(profx(i,:)*l_per_pix(1)*1000, 'k', 'linewidth', 1.2)
yline(profx(i,:)*l_per_pix(1)*1000, 'Color',cmapx(i,:))

end
cmapy = cool(size(profy,1))
for i = 1:size(profy,1)
xline(profy(i,:)*l_per_pix(2)*1000, 'k', 'linewidth', 1.2)
xline(profy(i,:)*l_per_pix(2)*1000, 'Color',cmapy(i,:))
end


t2 = nexttile(14, [1 5])
hold on
xSect= [];
for i = 1:size(profx,1)
    xSect(i,:) = mean(myheight(profx(i,1):profx(i,2),:), 'omitnan');
    plot(yPix, xSect(i,:),'Color',cmapx(i,:), 'linewidth', 1.5)
    plot(yPix, xSect(i,:), 'k')
end
for i = 1:length(locs)
    yline(locs(i), ':k')
end
ylim(heightLim)
xlim([min(yPix) max(xPix)])
xlabel('distance/µm')
ylabel('height/µm')

title('X section')
t2.Position(1) = t1.Position(1);
t2.Position(3) = t1.Position(3);


t3 = nexttile(1, [2 1])
hold on
ySect= [];
for i = 1:size(profy,1)
    ySect(i,:) = mean(myheight(:,profy(i,1):profy(i,2)), 2, 'omitnan');
    plot(ySect(i,:), xPix, 'Color', cmapy(i,:), 'linewidth', 1.5)
    plot(ySect(i,:), xPix, 'k')
end
for i = 1:length(locs)
    xline(locs(i), ':k')
end

axis tight
ylabel('distance/µm')
xlabel('height/µm')
ylim([min(xPix) max(xPix)])
xlim(heightLim)

title('X section')
linkaxes([t1 t2], 'x')
linkaxes([t1 t3], 'y')
colormap('hsv')


%% Height Histogram
% Flatten the height data and remove NaNs
heightVector = myheight(:);
heightVector = heightVector(~isnan(heightVector));

% Define bin edges in µm (adjust as needed)
binEdges = -3:1:10;  % Example: from -2 µm to 4 µm in 0.5 µm steps

% Compute histogram
[counts, edges] = histcounts(heightVector, binEdges);

% Convert counts to area
areaPercent = 100 * counts / sum(counts);

% Bin centers for plotting
binCenters = edges(1:end-1) + diff(edges)/2;

% Plot
figure;
bar(binCenters, areaPercent, 'FaceColor', [0.2 0.6 0.8]);
xlabel('Height Range (µm)');
ylabel('Area (%)');
title('Surface Area Percentage by Height Range');
grid on;

