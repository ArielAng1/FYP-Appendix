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

% NOTE: File location has to be on path set for MATLAB
% 
% myfile = uigetfile('*.*')
% fileID = fopen(myfile);
% A = fread(fileID,'ubit8');
% fclose(fileID);
% fileID = fopen(myfile);
% AA = fread(fileID,'ubit8');
% fclose(fileID);
% 
% Display name
% char(A(end-62:end)')
% 
% Fbase=[1 16^2 16^4 16^6]';
% for myioffset=[1:18]
%    i1=12+4*(myioffset-1)+1;
% myoffset(myioffset)=A(i1:i1+3)'*Fbase;
% end
% set equal to 1 for measurement conditions
% set equal to 2 for optical data
% set equal to 3 for laser optical data
% set equal to 4 for laser intensity data
% set equal to 7 for height data
% 
% Extract measurement conditions
% Extract length(pm) per pixel 
% typeoffset=1;
% MCbyte = [168 172 176];
% for i = 1:3
% l_per_pix(i) = typecast(uint8(A(myoffset(typeoffset)+MCbyte(i)+...
%     1:myoffset(typeoffset)+MCbyte(i)+4)), 'uint32');
% end
% l_per_pix = double(l_per_pix);
% 
% MCbyte = 148;
% 
% plComp = double(typecast(uint8(A(myoffset(typeoffset)+MCbyte+...
%     1:myoffset(typeoffset)+MCbyte+4)), 'uint32'));
% 
% Extract height data
% typeoffset=7;
% 
% next offset and reshape only the image matrix 
% 
% heightendoffset=min(myoffset((myoffset-myoffset(typeoffset))>0));
% C=A(myoffset(typeoffset)+797:heightendoffset)';
% heightrows=A(myoffset(typeoffset)+1:myoffset(typeoffset)+4)'*Fbase; 
% heightcols=A(myoffset(typeoffset)+5:myoffset(typeoffset)+8)'*Fbase; 
% heightcodingbase=A(myoffset(typeoffset)+9:myoffset(typeoffset)+12)'*Fbase; 
% largestheightallowed=A(myoffset(typeoffset)+29:myoffset(typeoffset)+32)'*Fbase;
% HF=reshape(C,4,heightrows*heightcols)'; 
% myheight=HF*Fbase;
% 
% output in micrometers
% myheight=flip(reshape(myheight,heightrows,heightcols)*l_per_pix(3)/1E6,1);
% myheight=myheight-mode(myheight(:));
% 
% construct pixel vectors in µm
% xPix = (1:size(myheight,1))*l_per_pix(1)/10^(3+plComp);
% yPix = (1:size(myheight,2))*l_per_pix(2)/10^(3+plComp);
% 
% Extract optical data
% 
% typeoffset=2;
% 
% opticalendoffset=min(myoffset((myoffset-myoffset(typeoffset))>0)); 
% C=A(myoffset(typeoffset)+21:opticalendoffset)'; 
% opticalrows=A(myoffset(typeoffset)+1:myoffset(typeoffset)+4)'*Fbase; 
% opticalcols=A(myoffset(typeoffset)+5:myoffset(typeoffset)+8)'*Fbase; 
% opticalcodingbase=A(myoffset(typeoffset)+9:myoffset(typeoffset)+12)'*Fbase; 
% largestopticalallowed=A(myoffset(typeoffset)+29:myoffset(typeoffset)+32)'*Fbase;
% OF=reshape(C,3,opticalrows*opticalcols)';
% myoptical=rescale(flip(flip(reshape(OF,[opticalrows opticalcols 3]),3),1));
% 
% Extract laser data
% typeoffset=3;
% 
% opticalendoffset=min(myoffset((myoffset-myoffset(typeoffset))>0)); 
% C=A(myoffset(typeoffset)+21:opticalendoffset)'; 
% opticalrows=A(myoffset(typeoffset)+1:myoffset(typeoffset)+4)'*Fbase; 
% opticalcols=A(myoffset(typeoffset)+5:myoffset(typeoffset)+8)'*Fbase; 
% opticalcodingbase=A(myoffset(typeoffset)+9:myoffset(typeoffset)+12)'*Fbase; 
% largestopticalallowed=A(myoffset(typeoffset)+29:myoffset(typeoffset)+32)'*Fbase;
% OF=reshape(C,3,opticalrows*opticalcols)';
% mylaser=rescale(flip(flip(reshape(OF,[opticalrows opticalcols 3]),3),1));

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


%% Nested grid search for tilt to give maximum peak in histogram

% mask
mask = 0;
xMask = [0 0.1; 0.9 1];
yMask = [0.0 0.1; 0.9 1];

xMask = round(xMask*xSize);
xMask(xMask==0) =1;
xMask(xMask>xSize) =xSize;

yMask = round(yMask*ySize);
yMask(yMask==0) =1;
yMask(yMask>ySize) =ySize;

% Set nested grid search parameters
gridsize = 3;
nestRatio = 3/2;
gridstep0 = 1E-1;
nestSize = 20;

% Initalise variables
slopex = [];
slopey = [];
peakH = [];
height = nan(size(myheight));
slopex(1) = 0;
slopey(1) = 0;
gridstep = gridstep0;

height = myheight+ wgn(size(myheight,1),size(myheight,2),-20);
%height = myheight(1:1:end, 1:1:end)+ wgn(length(myheight),1,-20);
tic
figure(20)
clf
for nest = 1:nestSize 
nest
gridstep = gridstep/nestRatio;
gridNest = (1:gridsize)*gridstep;
if nest == 1
gridx = gridNest -  mean(gridNest);
gridy = gridNest - mean(gridNest);
else
gridx = gridNest -  mean(gridNest)+slopex(nest-1);
gridy = gridNest - mean(gridNest)+slopey(nest-1);
end

for p = 1:gridsize
    for q = 1:gridsize
        if mask
            peakH(p,q) = peakHeightDist(gridx(p),gridy(q),...
                height([yMask(1):yMask(2), yMask(3):yMask(4)],[xMask(1):xMask(2), xMask(3):xMask(4)]));
        else
            peakH(p,q) = peakHeightDist(gridx(p),gridy(q),...
                height);
        end

    end
end

% nexttile
% image(peakH,'CDataMapping','scaled')

[p,q]=find(peakH==max(max(peakH)),1);

slopex(nest) = gridx(p);
slopey(nest) = gridy(q);
end
toc
nexttile
hold on
plot(slopex, 'o-')
plot(slopey, 'o-')
ylabel('slope (radians)')
xlabel('iternations')
legend(strcat("x tilt ", num2str(1000*gridx(p),3), " mrad"),...
    strcat("y tilt ", num2str(1000*gridy(q),3), " mrad")) 
nexttile
hold on
dslopex = abs(diff(slopex));
dslopey = abs(diff(slopey));
plot(dslopex, 'o-')
plot(dslopey, 'o-')
plot(gridstep0*gridsize/9./(nestRatio.^(1:nestSize)), 'o-')
set(gca, 'yscale', 'log')
ylabel('diff slope (radians)')
xlabel('iternations')
legend(strcat("x tilt ", num2str(1000*dslopex(end),3), " mrad"),...
    strcat("y tilt ", num2str(1000*dslopey(end),3), " mrad"),...
    strcat("tilt range ", num2str(1000*gridstep0*(gridsize-1)/2./(nestRatio^nestSize),3), " mrad")) 

sgtitle("Tilt iteration" )

% Apply tilt to data
myheight1 = myheight'+slopex(end)*(1:size(myheight,1));
myheight = myheight1'+slopey(end)*(1:size(myheight,2));

% set most common height as baseline
heightVector = sort(myheight(:));
heightVector =heightVector(~isnan(heightVector));
[f, xi] = ksdensity(heightVector);
myheight = myheight - xi(f == max(f));
heightLim = [heightVector(round(end/100)) heightVector(round(99*end/100))];

figure(5)
clf
hold on
h = surf(yPix, xPix, myheight,...
    'EdgeColor','none');



view([0 0 1]);

set(gca, 'Projection','perspective')
axis tight
xlabel('µm')
ylabel('µm')
zlabel('µm')
daspect([1 1 1]);
caxis(heightLim)
cb = colorbar
cb.Label.String = 'height/µm';

if mask
xline(xMask(:)*l_per_pix(1)*1000, 'r');
yline(yMask(:)*l_per_pix(2)*1000, 'r');
end

%% Display 3D data

% Vertical exaggeration
vertEx = 5; 
colormap('turbo')
figure(5)
% Vectorize the height data
heightVector = myheight(:);
% Remove NaNs
heightVector = heightVector(~isnan(heightVector));

% set limits to 1 and 99th percentile
heightSort = sort(heightVector);
heightLim = [heightSort(round(end/100)) heightSort(round(99*end/100))];


%   Height image 
h = surf(yPix, xPix, myheight,...
    'EdgeColor','none');
 
camlight('right');
camlight('left');
camlight('headlight');
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
vertEx = 50  ;

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
profx = round([0.10*ySize, 0.11*ySize]);

profy = round([0.81*xSize, 0.815*xSize]);


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


%%  Display full VK data
% Gaussian smooth
sz = size(myheight);

gaussK = 2;

% Crop ( = 0 for no cropping)
crop = 0;
xcrop = 1:800;
ycrop = 1:600;

%   profile in x/y direction: mean between pixel values in prof
%   xX = 1 for profile in x direction

% direction
xX = 0;

% bands
prof = round([0.15*sz(1), 0.25*sz(1);
    0.45*sz(2),0.55*sz(2)]);


% Vertical exaggeration
vertEx = 2;

figure(10)
clf

%   Height image 
nexttile
imshow(rescale(myheight'/max(max(myheight))));
colormap('jet')
view(-90,90)

for i = 1:size(prof,1)
    if xX
        yline(prof(i,1),'linewidth',i)
        yline(prof(i,2),'linewidth',i)
    else
        xline(prof(i,1),'linewidth',i)
        xline(prof(i,2),'linewidth',i)
    end
end

daspect([1 1 1]);
title('Height')

%   Histogram
nexttile
pd = fitdist(myheight(:),'kernel');
x = round(min(myheight(:))):0.1:round(max(myheight(:)));
y = pdf(pd,x);

plot(x,y)

xlim([-2 16])

myHeightPlot = flip(imgaussfilt(myheight,gaussK),2);

nexttile
hold on
xSect = [];
for i = 1:size(prof,1)
    if xX
        xSect(i,:) = mean(myheight(prof(i,1):prof(i,2),:));
        plot(xSect(i,:), 'k','linewidth',i)
    else
        xSect(i,:) = mean(myheight(:,prof(i,1):prof(i,2)),2); 
        plot(xSect(i,:), 'k','linewidth',i)
    end
end

axis tight
xlabel('µm')
ylabel('µm')

title('X section')
nexttile
hold on
for i = 1:size(prof,1)
plot((1:size(xSect,2))/size(xSect,2), sort(xSect(i,:),2), 'k','linewidth',i)
end
ylabel('height (µm)')

title('height cdf')

sgtitle({[myfile]},'Interpreter','none');



%%  Display
% Gaussian smooth
sz = size(myheight);

gaussK = 2;

% Crop ( = 0 for no cropping)
crop = 0;
xcrop = 1:800;
ycrop = 1:600;

%   profile in x/y direction: mean between pixel values in prof
%   xX = 1 for profile in x direction

% direction
xX = 0;

% bands
prof = round([0.4*sz(1), 0.6*sz(1);
    0.4*sz(2),0.6*sz(2)]);

% laser or optical for drape
laserDrape = 0;

% Vertical exaggeration
vertEx = 2;

% Brighten laser
brighten = 0;
if brighten > 0
    mylaser = imlocalbrighten(mylaser, brighten)
end

if crop
    mylaser = mylaser(xcrop, ycrop,:);
    myoptical = myoptical(xcrop, ycrop,:);
    myheight = myheight(xcrop, ycrop);
    xPix = xPix(xcrop);
    yPix = yPix(ycrop);
end
    

figure
%	Laser RGB image
subplot(2,7,1:2)
imshow(mylaser)
view(-90,90) 
title('Laser')

%	Optical RGB image
subplot(2,7,3:4)
imshow(myoptical)
view(-90,90) 
title('Optical')

%   Height image 
subplot(2,7,5:6)
imshow(rescale(myheight/max(max(myheight))));
colormap('jet')
view(-90,90)

for i = 1:size(prof,1)
    if xX
        yline(prof(i,1),'linewidth',i)
        yline(prof(i,2),'linewidth',i)
    else
        xline(prof(i,1),'linewidth',i)
        xline(prof(i,2),'linewidth',i)
    end
end

daspect([1 1 1]);
title('Height')

%   Histogram
subplot(2,7,7)
plot((1:length(myheight(:)))/length(myheight(:)), sort(myheight(:)), 'k');
ylabel('height (µm)')
axis tight
title('height cdf')

myHeightPlot = flip(imgaussfilt(myheight,gaussK),2);

% Drape image over height image 
subplot(2,7,8:11)
hold on
if laserDrape
    myDrape = imlocalbrighten(flip(mylaser,2));
    title(['3D Laser, vertical exaggeration: ', num2str(vertEx)])
else
    myDrape = flip(myoptical,2);
    title(strcat("3D Optical, vertical exaggeration: ", num2str(vertEx)))
end

h = surf(yPix, xPix, flip(imgaussfilt(myheight,2),2), 'FaceColor','texturemap',...
    'EdgeColor','none',...
    'Cdata',myDrape);
 
h.FaceLighting = 'gouraud';
% h.AmbientStrength = 0.8;
% h.DiffuseStrength = 0.8;
% h.SpecularStrength = 0.1;
% h.SpecularExponent = 25;


view([1 1 1]);
daspect([1 1 1/vertEx])
set(gca, 'Projection','perspective')
axis tight
xlabel('µm')
ylabel('µm')
zlabel('µm')

subplot(2,7,12:13)
hold on
xSect = [];
for i = 1:size(prof,1)
    if xX
        xSect(i,:) = mean(myheight(prof(i,1):prof(i,2),:));
        plot(yPix, xSect(i,:), 'k','linewidth',i)
    else
        xSect(i,:) = mean(myheight(:,prof(i,1):prof(i,2)),2); 
        plot(xPix, xSect(i,:), 'k','linewidth',i)
    end
end

axis tight
xlabel('µm')
ylabel('µm')

title('X section')
subplot(2,7,14)
hold on
for i = 1:size(prof,1)
plot((1:size(xSect,2))/size(xSect,2), sort(xSect(i,:),2), 'k','linewidth',i)
end
ylabel('height (µm)')

title('height cdf')

sgtitle({[myfile]},'Interpreter','none');

%% Add profile lines and zero level
z = get(gca, 'View');
az  = z(1);
if az>0&az<180
    col = 'red';
line(max(yPix)*ones(size(xPix)), xPix, zeros(size(xPix)), 'Color','black',...
    'LineWidth',0.5, 'LineStyle','-')
line(max(yPix)*ones(size(xPix)), xPix, myHeightPlot(:, end), 'Color', col,...
    'LineWidth', 2)
else
    line(min(yPix)*ones(size(xPix)), xPix, zeros(size(xPix)), 'Color','black',...
    'LineWidth',0.5, 'LineStyle','-')
line(min(yPix)*ones(size(xPix)), xPix, myHeightPlot(:, 1), 'Color', col,...
    'LineWidth', 2)
end

if az>-90&az<90
line(yPix, max(xPix)*ones(size(yPix)), zeros(size(yPix)), 'Color','black',...
    'LineWidth',0.5, 'LineStyle','-')
line(yPix, max(xPix)*ones(size(yPix)), myHeightPlot(1, :), 'Color', col,...
    'LineWidth', 2)

else
line(yPix, max(xPix)*ones(size(yPix)), zeros(size(yPix)), 'Color','black',...
    'LineWidth',0.5, 'LineStyle','-')
line(yPix, max(xPix)*ones(size(yPix)), myHeightPlot(end, :), 'Color', col,...
    'LineWidth', 2)
end

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


%%
function peakHeightDist = peakHeightDist(xSlope, ySlope, height)

myheight1 = height+ySlope*(1:size(height,2));
myheight1 = myheight1+xSlope*(1:size(height,1))';

xi = linspace(min(myheight1(:)), max(myheight1(:)), 2048);
y = myheight1(:);
y = sort(y(~isnan(y)));

% Kernel density estimation
%[f, ~] = ksdensity(y);
% Histogram
h = histogram(y);


peakHeightDist = max(h.Values);
end

