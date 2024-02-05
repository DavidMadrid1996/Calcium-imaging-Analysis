%% David Madrid . Last Rev 21/10/2023
clear
user_settings

%% Prompt initial user interface to load fish

% Checking the number of the trial 
temp_table(:,1) = table(list); % First column of the table --> list
temp_vector = 1:length(list); % Second column of the table --> 1:length(list)
temp_table(:,2) = table(temp_vector'); 
temp_table.Properties.VariableNames = {'List','nfish'};  
disp(temp_table) % Display the table as info
pause(0.5)

prompt = {'Enter nfish:','Enter sampling rate',... % 'Enter estimate CZI time' (2 Question)
    'Enter INITIAL Calcium baseframe(No activity)',...
    'Enter FINAL Calcium baseframe(No activity)','Enter apparatus','Enter region',...
    'Enter Indicator'};

dlgtitle = 'Check the table int eh Command Window';
dims = [1 35];

definput ={sprintf('%.0f',0),sprintf('%.0f',0),sprintf('%.0f',0),...
    sprintf('%.0f',0),'Confocal','Hindbrain','GCaMP6s'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

clear prompt

%% Loading fish parameters                    

z = str2double(answer(1)); % @SET number of the fish 
CALCIUM.estimate_time = str2double(answer(2)); % @SET
CALCIUM.indicator = char(answer(7)); % @ SET
CALCIUM.region= char(answer(6)); % @ SET
CALCIUM.apparatus= char(answer(5)); % @ SET
baseframe = str2double(answer(3)):str2double(answer(4)); % @SET! idx of frames to use as F0 in differential formula:
CALCIUM.baseframe = baseframe; 

%% Prompt initial user interface to z stack location

prompt = {'Most Dorsal ','Most Ventral ','(CURRENT) Location of the current slice in the Zstack'};

dlgtitle = 'Check the table in the Command Window';
dims = [1 35];

definput ={sprintf('%.0f',0),sprintf('%.0f',0),sprintf('%.0f',0)};

answer = inputdlg(prompt,dlgtitle,dims,definput);

DorsalMost = abs(str2double(answer(1))-str2double(answer(2))); 
temp = str2double(answer(3))-str2double(answer(2)); 
CurrentZ = (temp) / DorsalMost; 

CALCIUM.Zstack.Dorsal= str2double(answer(1));
CALCIUM.Zstack.Ventral = str2double(answer(2));
CALCIUM.Zstack.Current = str2double(answer(3));
CALCIUM.Zposition.DorsalRef = DorsalMost; 
CALCIUM.Zposition.CurrentRef= CurrentZ; 

clear temp clear answer CurrentZ DorsalMost

% 1 is when the Mth axons are crossing from the HB to the medulla. (DORSAL)
% 0 is when you see the Mth axons crossing. (VENTRAL)

%% Importing

% IMPORT a group list if there is one 
if isfile(fullfile(path.data,'grouplist.mat'))
    load (fullfile(path.data,'grouplist.mat'))  
else 
    
end 

nfish = z; % @ SET
CALCIUM.nfish = nfish ;

% Tiff_directory_loop =  fullfile(Tiff_directory , list(z)); % NOT IN USE
CSV_directory_loop = fullfile(CSV_directory , list2(z));
Video_directory_loop = fullfile(Videoreal_directory , list3(z));
CZI_directory_loop = fullfile(CZI_directory , list(z));
zstack_directory_loop = fullfile(zstack_directory  , 'zstack.czi');
% SNAP_directory_loop = fullfile(SNAP_directory  , list4(z));



% CALCIUM.Tiff_directory_loop = Tiff_directory_loop;  % NOT IN USE
CALCIUM.CSV_directory_loop = CSV_directory_loop ;
CALCIUM.Video_directory_loop = Video_directory_loop;
CALCIUM.CZI_directory_loop = CZI_directory_loop;
CALCIUM.zstack_directory_loop = zstack_directory_loop;
% CALCIUM.SNAP_directory_loop = SNAP_directory_loop; 

% Checking before if the video is good label 
CSV_file = table2array(readtable(CALCIUM.CSV_directory_loop,'NumHeaderLines',3)); %@ SET number of Headers
plot(CSV_file(:,15));
pause(0.5)

prompt = 'You want to continue with the analysis? y or n';
str = input(prompt,'s');
if isequal(str,'Y') || isequal(str,'y') 
elseif isequal(str,'N') || isequal(str,'n')
    return
end
close all


%% Number of file test in the csv folder,czi folder and video folder

CALCIUM.path = path;
CALCIUM.list = list;
CALCIUM.ref = char(CALCIUM.list(nfish,1));
grouplist{nfish} = strcat('CALCIUM_', char(CALCIUM.list(nfish,1))) ;
save(fullfile(path.data,'grouplist.mat'),'grouplist');
CALCIUM.trialref = CALCIUM.list(nfish,1); %save references for the trials
mkdir (fullfile(path.data,'dataCALCIUM',CALCIUM.ref)); % Creating the directory for dataCALCIUM
mkdir (fullfile(path.data,'datamovies',CALCIUM.ref)); % Creating the directory for datamovies
mkdir (fullfile(path.data,'datawaves',CALCIUM.ref));% Creating the directory for datawaves
mkdir (fullfile(path.data,'TIFF',CALCIUM.ref));
CALCIUMimg('save',CALCIUM,[],list,nfish);

%--------  If you have specific conditions
%{
for triali = 1:length(CALCIUM.trialref) % For now just one experimental condition
    CALCIUM.condition(triali,1) = C1(triali); 
    % CALCIUM.condition(triali,2) = CALCIUM.list(triali).C2;
    % CALCIUM.condition(triali,3) = CALCIUM.list(triali).C3;
    % ...
end

%}

%  CALCIUM.condition = C1(nfish);
CALCIUM.nonanidx=  1 ; %find(~isnan(CALCIUM.condition(:,1)));
 CALCIUMimg('save',CALCIUM,[],list,nfish);


% test saving
clearvars -except tiff nfish video csv list  
[CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);


                                                                      
%% 1.0 - IMPORT RAW MOVIE FROM IMAGES    

user_settings
data = bfopen(char(strcat(CALCIUM.CZI_directory_loop,'.czi')));
% czifinfo(char(strcat(CALCIUM.CZI_directory_loop,'.czi')))

for j = 1:length(data{1,1})
  movies4D (:,:,j) = data{1,1}{j,1};
end 

 movies4D = im2double(movies4D);
%  subplot(121)
%  imagesc(mean(movies4D,3)); colormap('jet')
 CALCIUM.mean = mean(movies4D,3) ;
%  subtitle('Brain Motion. Average trought time')
%  subplot(122)
%  imagesc(movies4D(:,:,1)); colormap('jet')
%  subtitle('Individual trial')
%  saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('Brain_Motion',char(CALCIUM.list(nfish,1)),'.fig')))
%  uiwait



% ------  If the automatic extracction doesn't work 
%{

% Save as sequence of images 

ds1 = imageDatastore(CALCIUM.Tiff_directory_loop); 
ds2 = imageDatastore(tiff);
%{
ds3 = imageDatastore('C:\Users\Usuario\OneDrive\Escritorio\211011\cerebellum_right_3\*.tif');
ds4 = imageDatastore('C:\Users\Usuario\OneDrive\Escritorio\4\*.tif');
%}

% Computing the matrix
sum = 0;
while hasdata(ds1) % If are single experiment you use just ds1 (3D matrix)

    
    img = read(ds1) ;             % read image from datastore
    % creates a new window for each image
    for i = 1
        sum = sum +1;
    end
    movies3D1(:,:,sum) = img; % changes in each while
end, movies4D_8(:,:,:,1) =  movies3D1; clear img movies3D1
sum = 0;
while hasdata(ds2)
    
    img = read(ds2) ;             % read image from datastore
    % creates a new window for each image
    for i = 1
        sum = sum +1;
    end
    movies3D2(:,:,sum) = img;
end, movies4D_8 (:,:,:,2) =  movies3D2; clear img movies3D2

%{ 
sum = 0;
while hasdata(ds3)
    
    img = read(ds3) ;             % read image from datastore
    % creates a new window for each image
    for i = 1
        sum = sum +1;
    end
    movies3D3(:,:,sum) = img;
end, movies4D_8 (:,:,:,3) =  movies3D3; clear img movies3D3
sum = 0;
while hasdata(ds4)
    
    img = read(ds4) ;             % read image from datastore
    % creates a new window for each image
    for i = 1
        sum = sum +1;
    end
    movies3D4(:,:,sum) = img;
end, movies4D_8 (:,:,:,4) =  movies3D4; clear img movies3D4
%}



%unit 8 to double
movies4D =  im2double(movies4D_8); clear movies4D_8
%}

movieref = '_00raw';
CALCIUMmov.ref = CALCIUM.ref;
CALCIUMmov.movieref= movieref;
CALCIUMmov.data = movies4D ;
CALCIUMmov.hist{1} = 'raw: imported + assemble4Dmovies';
% % % % ................... SAVING THE RAW MOVIE .......................% % % %
%CALCIUMimg('savemovie', CALCIUMmov, movieref,list,nfish);

% CALCIUM.srate = CALCIUM.estimate_time/size(movies4D,3); % @ SET can change everytime
CALCIUM.srate =CALCIUM.estimate_time; 
CALCIUM.timebase = CALCIUM.srate.*(1:size(movies4D,3));
CALCIUMimg('save',CALCIUM,[],list,nfish);  % If you want to save the raw movie

%% Importing zstack

if isempty(CALCIUM.zstack_directory_loop )
else
    data2 = bfopen(char(CALCIUM.zstack_directory_loop));

    for j = 1:length(data2{1,1})
        movies4D2 (:,:,j) = data2{1,1}{j,1};
    end
    movies4D2 = im2double(movies4D2);
    zstack = movies4D2; 
end

 CALCIUM.zstackImage = zstack(:,:,CALCIUM.Zstack.Current) ; 
% imagesc(histeq(CALCIUM.zstackImage)); colormap('bone')
 CALCIUMimg('save',CALCIUM,[],list,nfish);  % If you want to save the raw movie
 clear zstack

% 
% %% Importing snap
% 
% if isempty(CALCIUM.SNAP_directory_loop)
% else
%     data3 = bfopen(char(CALCIUM.SNAP_directory_loop));
% 
%     for j = 1
%         movies4D3 (:,:,j) = data3{1,1}{j,1};
%     end
%     movies4D3 = im2double(movies4D3);
%     SNAP = movies4D3; 
% 
% end
% 
% CALCIUM.SNAPImage = SNAP ; 
% CALCIUMimg('save',CALCIUM,[],list,nfish);  % If you want to save the raw movie
% clear SNAP
% 
% % Test
% %CALCIUMmov= CALCIUMimg('loadmovie',nfish,movieref,list,nfish);
% 
% % Checking sizes 
% % Carefull, apparently when the the aryscam simage is processed there is a
% % lost of the resolution. When a 1024x1024 image is processed the resulted
% % image is 1012x1012. To correct this, black pixels will be added in the
% % las columns and last rows. If the non processed image (correct resolution) 
% % is used instead of the processed one, the time of computation is two
% % high. Both, time series and SNAP must have the SAME resoltion to be able
% % to use the mask of the rois.
% 
% if size(CALCIUMmov.data(:,:,1)) == size(CALCIUM.SNAPImage)
% else 
%     error('Not the same resolution')
% end 
% 
% if size(CALCIUMmov.data(:,:,1)) == size(CALCIUM.SNAPImage)
% else 
% %     Blank = NaN(abs(size(CALCIUMmov.data(:,:,1),2)-size(CALCIUM.SNAPImage,2)),size(CALCIUM.SNAPImage,1)); 
% %     BlankRows = Blank; 
% %     BlankColumns  = Blank'; 
%     newMatrix = zeros(size(CALCIUM.SNAPImage,1),size(CALCIUM.SNAPImage,2),size(CALCIUMmov.data,3)); 
%     newMatrix(1:size(CALCIUMmov.data(:,:,1),1),1:size(CALCIUMmov.data(:,:,1),1),:) = CALCIUMmov.data; 
% clear CALCIUMmov.data
% CALCIUMmov.data = newMatrix; 
% size(CALCIUMmov.data); 
%     
% end 

 
%% NOT IN USE RIGHT NOW %%%%%%%%%%%%%%% BRAIN ALINEATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

%% 1.1 - BRAIN ALINEATION                                                                                              1.1
%    THROUGH TIME AND TRIALS ('_01registered')  


clearvars -except tiff nfish video csv list CALCIUM

% 1.1.1. REFERENCES for input/output movies (see 'z_notes.txt', point 5 for
% complete list)
inputRef =  '_00raw';
outputRef = '_01registered';
%load input movie (non-aligned raw)
[inputStruct] = CALCIUMimg('loadmovie',nfish,inputRef,list,nfish);
rawmov=inputStruct.data;

% 1.1.2. PERFORM COMPUTATIONS:  REGISTER
ref_frame = rawmov(:,:,1,CALCIUM.nonanidx(1)) ; %background from 1st 
% included trial
[registermov] = register_wrap(rawmov,ref_frame, CALCIUM.nonanidx); 
%{ 
% test
figure(2)
for i = 1:50
  imagesc( inputStruct.data(:,:,i,2))
  pause (0.1)
end 
for i = 1:50
  imagesc(registermov(:,:,i,2))
  pause (0.1)
end 
clf 
%}


% [Save reference frame in CALCIUM]
CALCIUM.backgr = ref_frame;
CALCIUMimg('save',CALCIUM,[],list,nfish);

% 1.1.3. SAVE NEW MOVIE STRUCTURE: Save new registered movie, copying some references from the movie
% structure used to apply new changes in
CALCIUMmov.ref = inputStruct.ref;
CALCIUMmov.movieref= outputRef;
CALCIUMmov.data = registermov;
CALCIUMmov.hist = inputStruct.hist;
CALCIUMmov.hist{length(CALCIUMmov.hist)+1,1} = 'register'; %append a new cell with new info
CALCIUMimg('savemovie', CALCIUMmov, CALCIUMmov.movieref,list,nfish);

clear inputStruct 

%}

%% NOT IN USE RIGHT NOW %%%%%%%%%%%%%%%%% DIFF WITHOUT PERCENTAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{ 

%% 1.2 - DIFFERENTIAL VALUES ('_02diff')                                                                     1.2
clearvars -except CALCIUM nfish

% 1.2.1. REFERENCES for input/output movies
inputRef = '_00raw'; %@ If you are using Brain Aligniation change '_00raw' for '_01registered'
outputRef = '_02diff';

[inputStruct] = CALCIUMimg('loadmovie',nfish,inputRef);

% 1.2.2. PERFORM COMPUTATIONS: %DIFFERENTIAL VALUES

inputdata = inputStruct.data;

baseframe = 1; % @SET! idx of frames to use as F0 in differential formula
% Turn into string to save later in History:
baseltext = strcat(num2str(baseframe(1)),'to',num2str(baseframe(end)));

% Preallocate in NaN
inputdim = size(inputdata);
diffmovies = NaN(inputdim(1),inputdim(2),inputdim(3)+1,inputdim(4));

for triali = makeRow(CALCIUM.nonanidx) %import only included trials
    inputmovie = squeeze(inputdata(:,:,:,triali));
    diffmovies(:,:,:,triali) = raw2diff(inputmovie, baseframe);
    CALCIUM.backgr(:,:,triali) = diffmovies(:,:,end,triali); % store background
end

% 1.2.3. SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
% structure used to apply new changes in
CALCIUMmov.ref = inputStruct.ref;
CALCIUMmov.movieref= outputRef;
CALCIUMmov.data = diffmovies;
CALCIUMmov.hist = inputStruct.hist;
CALCIUMmov.hist{length(CALCIUMmov.hist)+1,1} = outputRef; 
%append a new cell with new info

CALCIUMimg('savemovie', CALCIUMmov, CALCIUMmov.movieref);



% SUGGESTION: if different F0 are ,
%keep the basic reference + info about the F0,
% e.g. outputRef = '_02diffbase10';


%}


%% 1.2 - DIFFERENTIAL VALUES PERCENTAGE ('_03diff_perc')                                 1.2

clearvars -except tiff nfish video csv list CALCIUM CALCIUMmov

% For autoloading after saving
%{

% 1. REFERENCES for input/output movies
inputRef =  '_00raw'; %@ If you star using more than one trial use '_01registered'
outputRef = '_03diff_perc';

[inputStruct] = CALCIUMimg('loadmovie',nfish,inputRef,list,nfish);

% 2. PERFORM COMPUTATIONS: %DIFFERENTIAL VALUES

inputdata = inputStruct.data;

%}

inputStruct =  CALCIUMmov; 
inputdata = CALCIUMmov.data;
outputRef = '_03diff_perc';


% Turn into string to save later in History:
baseltext = strcat(num2str(CALCIUM.baseframe(1)),'to',num2str(CALCIUM.baseframe(end)));


% Preallocate in NaN
inputdim = size(inputdata);
%diffmovies = NaN(inputdim(1),inputdim(2),inputdim(3)+1,inputdim(4)); for 4D
diffmovies = NaN(inputdim(1),inputdim(2),inputdim(3)+1); % for 3D

%For 4D
%{  

for triali = makeRow(CALCIUM.nonanidx) %import only included trials
    inputmovie = squeeze(inputdata(:,:,:,triali));
    diffmovies(:,:,:,triali) = raw2diffperc2(inputmovie, baseframe);
    CALCIUM.backgr(:,:,triali) = diffmovies(:,:,end,triali); % store background
end

%} 

inputmovie = inputdata(:,:,:);
diffmovies(:,:,:) = raw2diffperc2(inputmovie, CALCIUM.baseframe);
CALCIUM.backgr(:,:,:) = diffmovies(:,:,end); % store background

% 3.SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
% structure used to apply new changes in
CALCIUMmov.ref = inputStruct.ref;
CALCIUMmov.movieref= outputRef;
CALCIUMmov.data = diffmovies;
CALCIUMmov.hist = inputStruct.hist;
CALCIUMmov.hist{length(CALCIUMmov.hist)+1,1} = outputRef; %append a new cell with new info

% % % % ................... SAVING THE DIFF MOVIE .......................% % % %
%CALCIUMimg('savemovie', CALCIUMmov, CALCIUMmov.movieref,list,nfish);


% SUGGESTION: if different F0 are ,
%keep the basic reference + info about the F0,
% e.g. outputRef = '_02diffbase10';


%{
%% 1.3- [1] MAKE CROPMASK,                                                                                          1.3
%     [2] CROPPED-MOVIES (optional)                           

clearvars -except tiff nfish video csv list CALCIUM sr 

% [1] MAKE & SAVE CROP MASK

% MAKE MASK FROM REFERENCE FRAME AND SAVE IN CALCIUM
ref_frame = CALCIUM.backgr(:,:,CALCIUM.nonanidx(1)); 
%the background from the first included trial

%  Before cropping: check all backgrounds to take into account if there is
%  much movements (and, for instance, leave out of the mask the margins)
%{
for triali = makeRow(CALCIUM.nonanidx)
imagesc(CALCIUM.backgr(:,:,triali)); 
title(strcat('trial=',num2str(triali)))
    pause 
    %to advance to the next frame, 
    %press any key; to skip to the end, press 'Ctrl+C'
end
%}

% DRAW & SAVE CROPMASK:
[crop_poly, crop_mask] = roi_draw(ref_frame);
close 

% View the result on a selected trial
trialsel = 1; %@ SET (if you want to check the mask onto any specific frame)
roi_preview(CALCIUM.backgr(:,:,trialsel), crop_poly{1}); 
saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('Global_signal_',char(CALCIUM.list(nfish,1)))))
close 
% polygon_preview(CALCIUM.backgr(:,:,1), CALCIUM.crop.poly{1,1}); 

% IF WE ARE HAPPY WITH THE MASK: SAVE in structure
CALCIUM.crop.mask = crop_mask; 
CALCIUM.crop.poly = crop_poly{1}; %stored in rows

% Save one cropped image to help in drawing the different ROI
cropframe =roi_crop(CALCIUM.backgr(:,:,CALCIUM.nonanidx(1)), CALCIUM.crop.mask);

CALCIUMimg('save',CALCIUM,[],list,nfish);

% 
% [2] CROPPED-MOVIES (mute if you don't want to extract them)

% 1. REFERENCES for input/output movies
inputRef = '_03diff_perc';
% inputRef =  '_03diff_perc';
outputRef = '_04crop';

inputStruct = CALCIUMimg('loadmovie',nfish,inputRef,list,nfish);
inputdata = inputStruct.data;

% 2. PERFORM COMPUTATIONS: APPLY MASK
cropmovies = NaN (size(inputdata));

%initialize
for triali = makeRow(CALCIUM.nonanidx)
inputmovie = squeeze(inputdata(:,:,:,triali));
cropmovies(:,:,:,triali)= roi_crop(inputmovie, CALCIUM.crop.mask);
end
CALCIUM.crop.preview = cropmovies(:,:,end,CALCIUM.nonanidx(1)); 
CALCIUMimg('save',CALCIUM,[],list,nfish);

% 3.SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
% structure used to apply new changes in
CALCIUMmov.ref = inputStruct.ref;
CALCIUMmov.movieref= outputRef;
CALCIUMmov.data = cropmovies;  
CALCIUMmov.hist = inputStruct.hist;
CALCIUMmov.hist{length(CALCIUMmov.hist)+1,1} = 'crop_backgr'; %append a new cell with new info

CALCIUMimg('savemovie', CALCIUMmov, CALCIUMmov.movieref,list,nfish); 


%}

%{

%% s02 FILTERING                                                                      %%%%%%%%%%%%%%
%  (further preprocessing) and crop                                     %%%%%%%%%%%%%%
                                                                                                         %%%%%%%%%%%%%%
clearvars -except tiff nfish video csv list CALCIUM sr CALCIUMmov
user_settings
[CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);
% FILT1:

%@ SET FILTERS PARAMETERS :
% tcnst = 10;% Time-constant for TEMPORAL SMOOTHING
mean_filter = 9; %  gaussian smoothing kernel for SPATIAL SMOOTHERING
% meanpix = 3;
 medianpix = 9; % pixels size for MEDIAN SPATIAL FILTER


% 1. REFERENCES for input/output movies (see 'z_notes.txt', point 5 for
% complete list)
inputRef =   '_03diff_perc';
outputRef = '_05filt1'; %@ SET

% Load input movie
% For autoloading after saving
% [inputStruct] = CALCIUMimg('loadmovie',nfish,inputRef,list,nfish);

[inputStruct]=CALCIUMmov; 
inputdata=inputStruct.data;


% 2. PERFORM COMPUTATIONS: %DIFFERENTIAL VALUES

% Preallocate in NaN
filtmov = NaN(size(inputdata));

% % 2.1. Tcnst % Don't apply !
% for triali = makeRow(CALCIUM.nonanidx)
%     tempmov = inputdata(:,:,:,triali);
%     filt1(:,:,:,triali) = filter_Tcnst(tempmov,tcnst);
%     clear tempmov
% end


% 2.2. Spatial Filter (mean)
for triali = makeRow(CALCIUM.nonanidx)
    tempmov = inputdata (:,:,:,triali);
    filt2(:,:,:,triali) = filter_spatial(tempmov, mean_filter);
    clear tempmov
end

% 2.3. Median spatial filter
for triali = makeRow(CALCIUM.nonanidx)
    tempmov = filt2(:,:,:,triali);
    filt3(:,:,:,triali) = filter_median(tempmov, medianpix);
    clear tempmov
end

% % 2.4. Cubic filter
% for triali = makeRow(CALCIUM.nonanidx)
%     tempmov = filt3(:,:,:,triali);
%     filt4(:,:,:,triali) = filter_cubicBV(tempmov);
%     clear tempmov
% end

% % 2.5. Crop background
% for triali = makeRow(CALCIUM.nonanidx)
%     tempmov = filt4(:,:,:,triali);
%     filtmov(:,:,:,triali)= roi_crop(tempmov, CALCIUM.crop.mask);
%     clear tempmov
% end

filtmov = filt3; % we save the last filtered movie

% 3.SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
% structure used to apply new changes in
CALCIUMmov.ref = inputStruct.ref;
CALCIUMmov.movieref= outputRef;
CALCIUMmov.data = filtmov;

%@ SET !!! according to the filters applie (append as many as needeed)
CALCIUMmov.hist = inputStruct.hist;
%CALCIUMmov.hist{length(CALCIUMmov.hist)+1,1} = 'tcnst = 10'; %@ SET
CALCIUMmov.hist{length(CALCIUMmov.hist)+1,1} = 'spatialmean = 3'; %@ SET
CALCIUMmov.hist{length(CALCIUMmov.hist)+1,1} = 'median = 3'; %@ SET
% CALCIUMmov.hist{length(CALCIUMmov.hist)+1,1} = 'cubicBV'; %@ SET
% CALCIUMmov.hist{length(CALCIUMmov.hist)+1,1} = 'crop-background'; %@ SET
% CALCIUMimg('savemovie', CALCIUMmov, CALCIUMmov.movieref,list,nfish); 



%}

%% s03 DEFINE ROIS                                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                               
% Draw ROIs from neurons.
clearvars -except tiff nfish video csv list  CALCIUMmov
close all
user_settings
[CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);
 
% Preview to count neurons. Click twice in the last one.
figure('units','normalized','outerposition',[0 0 1 1])
 imagesc(CALCIUM.mean), colormap('jet'); % Use in you dont have the SNAP
% imagesc(adapthisteq(CALCIUM.SNAPImage)), colormap('bone'); % Use if you have the SNAP
disp('Mark the neurons that you want to analyze. Double click in the last one.')
clear num_neurons_x num_neurons_y;[num_neurons_x num_neurons_y] = getpts;
close all
figure('units','normalized','outerposition',[0 0 1 1])
imagesc(CALCIUM.mean);colormap ('jet'), hold on
% imagesc(adapthisteq(CALCIUM.SNAPImage)), colormap('bone'); hold on
scatter(num_neurons_x,num_neurons_y)
title('Labelled neurons')
saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('Labelled_neurons_',char(CALCIUM.list(nfish,1)),'.fig')))
close all
X = ['Number of neurons selected:',num2str(length(num_neurons_x))];
disp(X)


% Dynamic ROIs names creation and store in structure. The number of ROIs is base in the previous varibale
% length(num_neurons); 'roi_labels': in the order in which they will be plot (}%). Exmaple:
for i = 1:length(num_neurons_x)
    CALCIUM.roi.labels {i} = strcat('N',num2str(i)); 
end 
display(CALCIUM.roi.labels)
close all



% DEFINE with the help of the function. You can zoom in the image for a
% better ROI selection
% [manual_poly, manual_mask] = roi_draw(CALCIUM.mean,CALCIUM.roi.labels,num_neurons_x,num_neurons_y); %in each column, one set of rois
[manual_poly, manual_mask] = roi_draw2(CALCIUM.mean,CALCIUM.roi.labels,num_neurons_x,num_neurons_y,CALCIUMmov.data); %in each column, one set of rois
close all
% View the result: 1. Plot just the ROIs; 2. Plot ROIs + Names = {N1, N2, ...}
roi_preview_multiple(CALCIUM.mean, manual_poly, CALCIUM);
% roi_preview_multiple(adapthisteq(CALCIUM.SNAPImage), manual_poly, CALCIUM);


% C = imfuse(CALCIUM.mean,CALCIUM.SNAPImage,'falsecolor');
% imshow(C)

saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('ROIS_',char(CALCIUM.list(nfish,1)),'.fig')))
uiwait

close all

% Saving the mask and coordinates in the CALCIUM structure
CALCIUM.num_neurons_x = num_neurons_x;
CALCIUM.num_neurons_y = num_neurons_y;
CALCIUM.neurons = length(num_neurons_x) ; % Number of selected neurons @ SET
CALCIUM.roi.manual_poly  = manual_poly; % Coordinates
CALCIUM.roi.manual_mask = manual_mask; % Mask
CALCIUMimg('save',CALCIUM,[],list,nfish);

clear num_neurons_x num_neurons_y manual_poly manual_mask

%%%% Mid reference. Just one elongated ROI in the Midline

pause(0.2)

prompt = {'Do you want to have a medial reference?'};
dlgtitle = 'Input';
dims = [1 35];

definput ={'Y'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
answer = char(answer); 

if isequal(answer,'Y')
    figure('units','normalized','outerposition',[0 0 1 1])
    imagesc(CALCIUM.mean), colormap('jet');
%     imagesc(adapthisteq(CALCIUM.SNAPImage)), colormap('bone');
    title('Mark the mid reference/s')
    disp('Mark the mid reference/s.')
    clear num_neurons_x num_neurons_y;[num_neurons_x num_neurons_y] = getpts;
    close all
    imagesc(CALCIUM.mean);colormap ('jet'), hold on
%     imagesc(adapthisteq(CALCIUM.SNAPImage)), colormap('bone'); hold on
    scatter(num_neurons_x,num_neurons_y)
    title('Mid reference/s')
    saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('Mid_reference_',char(CALCIUM.list(nfish,1)),'.fig')))
    close all
    X = ['Number of mid reference/s:',num2str(length(num_neurons_x))];
    disp(X)

%Dynamic ROIs names creation and store in structure. The number of ROIs is base in the previous varibale
% length(num_neurons); 'roi_labels': in the order in which they will be plot (}%). Exmaple:
for i = 1:length(num_neurons_x)

    CALCIUM.Midreference.labels {i} = strcat('MID',num2str(i)); 
end 
display(CALCIUM.Midreference.labels)
close all



% DEFINE with the help of the function. You can zoom in the image for a
% better ROI selection
[manual_poly, manual_mask] = roi_draw(CALCIUM.mean,CALCIUM.Midreference.labels,num_neurons_x,num_neurons_y);
% [manual_poly, manual_mask] = roi_draw(adapthisteq(CALCIUM.SNAPImage),CALCIUM.Midreference.labels,...
%     num_neurons_x,num_neurons_y); %in each column, one set of rois
close all

% View the result: 1. Plot just the ROIs; 2. Plot ROIs + Names = {N1, N2, ...}
roi_preview_multiple(CALCIUM.mean, manual_poly, CALCIUM);
% roi_preview_multiple(adapthisteq(CALCIUM.SNAPImage), manual_poly, CALCIUM);
saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('Mid_reference',char(CALCIUM.list(nfish,1)),'.fig')))
close all

% Saving the mask and coordinates in the CALCIUM structure
CALCIUM.num_neurons_x_MID = num_neurons_x;
CALCIUM.num_neurons_y_MID = num_neurons_y;
CALCIUM.midreference= length(num_neurons_x) ; % Number of selected neurons @ SET
CALCIUM.Midreference.manual_poly  = manual_poly; % Coordinates
CALCIUM.Midreference.manual_mask = manual_mask; % Mask
CALCIUMimg('save',CALCIUM,[],list,nfish);

clear num_neurons_x num_neurons_y manual_poly manual_mask

elseif isequal(answer,'N')
else 
    error('Answer must be Y or N');

end


%%%% Lateral reference
% 
% pause(0.2)
% 
% prompt = {'Do you want to have a lateral reference?'};
% dlgtitle = 'Input';
% dims = [1 35];
% 
% definput ={'Y'};
% answer = inputdlg(prompt,dlgtitle,dims,definput);
% answer = char(answer); 
% 
% if isequal(answer,'Y')
%     figure('units','normalized','outerposition',[0 0 1 1])
%     imagesc(CALCIUM.mean), colormap('jet');
% %     imagesc(histeq(CALCIUM.SNAPImage)), colormap('bone')
%     title('Mark the lateral reference/s')
%     disp('Mark the lateralreference/s.')
%     clear num_neurons_x num_neurons_y;[num_neurons_x num_neurons_y] = getpts;
%     close all
%      imagesc(CALCIUM.mean);colormap ('jet'), hold on
% %     imagesc(histeq(CALCIUM.SNAPImage)), colormap('bone'), hold on
%     scatter(num_neurons_x,num_neurons_y)
%     title('lateral reference/s')
%     saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('Mid_reference_',char(CALCIUM.list(nfish,1)),'.fig')))
%     close all
%     X = ['Number of lateral reference/s:',num2str(length(num_neurons_x))];
%     disp(X)
% 
% %Dynamic ROIs names creation and store in structure. The number of ROIs is base in the previous varibale
% % length(num_neurons); 'roi_labels': in the order in which they will be plot (}%). Exmaple:
% for i = 1:length(num_neurons_x)
% 
%     CALCIUM.Latreference.labels {i} = strcat('LAT',num2str(i)); 
% end 
% display(CALCIUM.Latreference.labels)
% close all
% 
% 
% 
% % DEFINE with the help of the function. You can zoom in the image for a
% % better ROI selection
% [manual_poly, manual_mask] = roi_draw(CALCIUM.mean,CALCIUM.Latreference.labels,num_neurons_x,num_neurons_y);
% % [manual_poly, manual_mask] = roi_draw(histeq(CALCIUM.SNAPImage),CALCIUM.Latreference.labels,...
% %     num_neurons_x,num_neurons_y); %in each column, one set of rois
% close all
% 
% % View the result: 1. Plot just the ROIs; 2. Plot ROIs + Names = {N1, N2, ...}
% roi_preview_multiple(CALCIUM.mean, manual_poly, CALCIUM);
% % roi_preview_multiple(histeq(CALCIUM.SNAPImage), manual_poly, CALCIUM);
% saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('lateral_reference',char(CALCIUM.list(nfish,1)),'.fig')))
% 
% % Saving the mask and coordinates in the CALCIUM structure
% CALCIUM.num_neurons_x_LAT = num_neurons_x;
% CALCIUM.num_neurons_y_LAT = num_neurons_y;
% CALCIUM.latreference= length(num_neurons_x) ; % Number of selected neurons @ SET
% CALCIUM.Latreference.manual_poly  = manual_poly; % Coordinates
% CALCIUM.Latreference.manual_mask = manual_mask; % Mask
% CALCIUMimg('save',CALCIUM,[],list,nfish);
% 
% clear num_neurons_x num_neurons_y manual_poly manual_mask
% 
% elseif isequal(answer,'N')
% else 
%     error('Answer must be Y or N');
% 
% end
% close all 


%%%% References for the zstack

%%%% Lateral reference

pause(0.2)

prompt = {'Do you want to have a zstack lateral reference?'};
dlgtitle = 'Input';
dims = [1 35];

definput ={'Y'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
answer = char(answer); 

if isequal(answer,'Y')
    figure('units','normalized','outerposition',[0 0 1 1])
%   imagesc(CALCIUM.mean), colormap('jet');
    imagesc(histeq(CALCIUM.zstackImage)), colormap('bone'), hold on
    title('Mark the lateral reference/s')
    disp('Mark the lateralreference/s.')
    clear num_neurons_x num_neurons_y;[num_neurons_x num_neurons_y] = getpts;
    close all
%     imagesc(CALCIUM.mean);colormap ('jet'), hold on
imagesc(histeq(CALCIUM.zstackImage)), colormap('bone'), hold on
    scatter(num_neurons_x,num_neurons_y)
    title('lateral reference/s')
    saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('zstacklateral_reference_',char(CALCIUM.list(nfish,1)),'.fig')))
    close all
    X = ['Number of lateral reference/s:',num2str(length(num_neurons_x))];
    disp(X)

%Dynamic ROIs names creation and store in structure. The number of ROIs is base in the previous varibale
% length(num_neurons); 'roi_labels': in the order in which they will be plot (}%). Exmaple:
for i = 1:length(num_neurons_x)

    CALCIUM.ZstackLatreference.labels {i} = strcat('LAT',num2str(i)); 
end 
display( CALCIUM.ZstackLatreference.labels)
close all



% DEFINE with the help of the function. You can zoom in the image for a
% better ROI selection
[manual_poly, manual_mask] = roi_draw(histeq(CALCIUM.zstackImage),CALCIUM.ZstackLatreference.labels,...
    num_neurons_x,num_neurons_y); %in each column, one set of rois
close all

% View the result: 1. Plot just the ROIs; 2. Plot ROIs + Names = {N1, N2, ...}
% roi_preview_multiple(CALCIUM.mean, manual_poly, CALCIUM);
roi_preview_multiple(histeq(CALCIUM.zstackImage), manual_poly, CALCIUM);
saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('zstacklateral_reference',char(CALCIUM.list(nfish,1)),'.fig')))

% Saving the mask and coordinates in the CALCIUM structure
CALCIUM.num_neurons_x_ZLAT = num_neurons_x;
CALCIUM.num_neurons_y_ZLAT = num_neurons_y;
CALCIUM.Zlatreference= length(num_neurons_x) ; % Number of selected neurons @ SET
CALCIUM.ZLatreference.manual_poly  = manual_poly; % Coordinates
CALCIUM.ZLatreference.manual_mask = manual_mask; % Mask
CALCIUMimg('save',CALCIUM,[],list,nfish);

clear num_neurons_x num_neurons_y manual_poly manual_mask

elseif isequal(answer,'N')
else 
    error('Answer must be Y or N');

end

%%%% zstack medial reference

pause(0.2)

prompt = {'Do you want to have a zstack medial reference?'};
dlgtitle = 'Input';
dims = [1 35];

definput ={'Y'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
answer = char(answer); 

if isequal(answer,'Y')
    figure('units','normalized','outerposition',[0 0 1 1])
%   imagesc(CALCIUM.mean), colormap('jet');
    imagesc(histeq(CALCIUM.zstackImage)), colormap('bone'),hold on;
    title('Mark the medial reference/s')
    disp('Mark the medialreference/s.')
    clear num_neurons_x num_neurons_y;[num_neurons_x num_neurons_y] = getpts;
    close all
%     imagesc(CALCIUM.mean);colormap ('jet'), hold on
imagesc(histeq(CALCIUM.zstackImage)), colormap('bone'),hold on;
    scatter(num_neurons_x,num_neurons_y)
    title('medial reference/s')
    saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('zstackMid_reference_',char(CALCIUM.list(nfish,1)),'.fig')))
    close all
    X = ['Number of medial reference/s:',num2str(length(num_neurons_x))];
    disp(X)

%Dynamic ROIs names creation and store in structure. The number of ROIs is base in the previous varibale
% length(num_neurons); 'roi_labels': in the order in which they will be plot (}%). Exmaple:
for i = 1:length(num_neurons_x)

    CALCIUM.ZstackMidreference.labels {i} = strcat('LAT',num2str(i)); 
end 
display( CALCIUM.ZstackMidreference.labels)
close all



% DEFINE with the help of the function. You can zoom in the image for a
% better ROI selection
[manual_poly, manual_mask] = roi_draw(histeq(CALCIUM.zstackImage),CALCIUM.ZstackMidreference.labels,...
    num_neurons_x,num_neurons_y); %in each column, one set of rois
close all

% View the result: 1. Plot just the ROIs; 2. Plot ROIs + Names = {N1, N2, ...}
% roi_preview_multiple(CALCIUM.mean, manual_poly, CALCIUM);
roi_preview_multiple(histeq(CALCIUM.zstackImage), manual_poly, CALCIUM);
saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('zstackMid_reference',char(CALCIUM.list(nfish,1)),'.fig')))

% Saving the mask and coordinates in the CALCIUM structure
CALCIUM.num_neurons_x_ZMID = num_neurons_x;
CALCIUM.num_neurons_y_ZMID = num_neurons_y;
CALCIUM.ZMidreference= length(num_neurons_x) ; % Number of selected neurons @ SET
CALCIUM.ZMIDreference.manual_poly  = manual_poly; % Coordinates
CALCIUM.ZMIDreference.manual_mask = manual_mask; % Mask
CALCIUMimg('save',CALCIUM,[],list,nfish);

clear num_neurons_x num_neurons_y manual_poly manual_mask

elseif isequal(answer,'N')
else 
    error('Answer must be Y or N');

end
close all 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% CREATE WAVES STRUCRURE (for extraction in '_extract_LOItimeseries):

CALCIUMroiTS.ref =CALCIUM.ref; 
CALCIUMroiTS.roi =CALCIUM.roi; 
CALCIUMimg('savewave', CALCIUMroiTS,[],list,nfish);

%% s04 EXTRACT ROI timeseries                               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except tiff nfish video csv list CALCIUMmov
close
user_settings
[CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);

CALCIUMroiTS =CALCIUMimg('loadwave',nfish,[],list,nfish);

% CREATE STRUCTURES FROM SELECTED MOVIES

% 1. REFERENCES for input movies/output waves (see 'z_notes.txt', point 5 for
% complete list)
inputRef =  '_03diff_perc'; %@ If you star using filtered use '_05filt1'
fieldref = strcat(inputRef(4:end),inputRef(2:3)); 

%load input movie 
%For autoloading
% [inputStruct] = CALCIUMimg('loadmovie',nfish,inputRef,list,nfish);
[inputStruct] = CALCIUMmov;
movies=inputStruct.data;

% 2. PERFORM COMPUTATIONS:  EXTRACT WAVES 
nroi =length(CALCIUMroiTS.roi.labels);
allroi_waves = NaN(length(CALCIUM.timebase), nroi, CALCIUM.nonanidx(end)); 

for triali = makeRow(CALCIUM.nonanidx)
    movietrial = movies(:,:,:,triali);
    for roii = 1:nroi
        roimask = CALCIUMroiTS.roi.manual_mask(:,:,roii);
        allroi_waves(:,roii,triali) = roi_TSave(movietrial,roimask);
    end
end

% 3.SAVE in TS STRUCTURE: Save waves (CALCIUMroiTS) structure copying some references from the movie
% structure used to apply new changes in
CALCIUMroiTS.(fieldref).data= allroi_waves;
CALCIUMroiTS.(fieldref).times = CALCIUM.timebase;
CALCIUMroiTS.(fieldref).hist = inputStruct.hist;
CALCIUMimg('savewave', CALCIUMroiTS, CALCIUMroiTS.ref,list,nfish);
clear inputStruct

%% s05 VISUALIZE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except tiff nfish video csv list 
user_settings
[CALCIUM] = CALCIUMimg('load',nfish,[],list,nfish);

% 1. WAVEPLOTS OF TRIALS for each roi
% Load data
CALCIUMroiTS = CALCIUMimg('loadwave',nfish,[],list,nfish);
waves = CALCIUMroiTS.diff_perc03.data; %@ SET
% 
% % % Make waveplots
% % title = strcat('fish',num2str(CALCIUMroiTS.ref)); %@ SET
% % close all
% % plot_wavesplot(waves, CALCIUMroiTS.roi.labels, title ,CALCIUM.nonanidx, 2);
% 
% % 2. ERPs
% 
% %@ SET: what to plot:
% timeseries= CALCIUMroiTS.diff_perc03;
% roi2plotidx = 1:length(CALCIUMroiTS.roi.labels);
% titletxt = strcat('#',num2str(CALCIUMroiTS.ref));
% trials2plot =find(CALCIUM.condition==1); %@ SET 
% % trials2plot = choosetrials(CALCIUMI.condition, [1 2]);
% stim_dur = 0.0; %ms
% 
% %plot_ERPs(-timeseries,CALCIUM.timebase,trials2plot, roi2plotidx , CALCIUMroiTS.roi.labels, titletxt, stim_dur); 
% % Subplot_ERPs(-timeseries,CALCIUM.timebase,trials2plot, roi2plotidx , CALCIUMroiTS.roi.labels, titletxt, stim_dur,'ROI WAVE PLOTS'); 
% 
% % saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('waves_',char(CALCIUM.list(nfish,1)))))
% % close all


%% DEEPLABCUT ANALYSIS IMPORT

% CSV directory
CSV_file = table2array(readtable(CALCIUM.CSV_directory_loop,'NumHeaderLines',3)); %@ SET number of Headers
% Frame vector
frames = CSV_file(1:end-1,1); % @ SET sampling rate. The last element is removed
% because in the excel file the coordinates are 0 for that frame, so is not
% included apparently



%% Laser detection
light.likehood = round(CSV_file(:,4));
laserframeON = find( light.likehood,1,'first');
laserframeOFF = find( light.likehood,1,'last');

laserframe = laserframeON:laserframeOFF;


% Test 
% Read Video Object
vidObj = VideoReader(CALCIUM.Video_directory_loop);
time = (1/vidObj.FrameRate).*(frames+1); % Seconds. This is the complete duration of the video
lasertime = time(laserframeON:laserframeOFF); % This is the time vector when there is 
% laser light in the confocal and is being recorded by the camera
CALCIUMroiTS.deeplabcut.sr = vidObj.FrameRate;
% Create 4D a matrix with frames ON
framesON = read(vidObj,[1 laserframeON+1]);

% Use just one channel and go from 4D to 3D
frames_laserON = squeeze(framesON(:,:,2,:));
% Checking framer before, during and after laser
 subplot(231)
imagesc(frames_laserON(:,:,laserframeON-1))
title('Bfr laserON')

 subplot(232)
imagesc(frames_laserON(:,:,laserframeON))
title('laserframeON')

 subplot(233)
imagesc(frames_laserON(:,:,laserframeON+1))
title('Aft laserON')







saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('Laser_test_',char(CALCIUM.list(nfish,1)),'.fig')))

uiwait

%% Body Parts. Must be modified if there is a change in the config. project of deeplabcut

time_sync = time(1:length(lasertime));

CALCIUMroiTS.deeplabcut.light.likehood= light.likehood;
CALCIUMroiTS.deeplabcut.laserframe = laserframe; 
CALCIUMroiTS.deeplabcut.lasertime = lasertime;
CALCIUMroiTS.deeplabcut.laserframeON  = laserframeON ;
CALCIUMroiTS.deeplabcut.laserframeOFF  = laserframeOFF ;
CALCIUMroiTS.deeplabcut.lasertimesync = time_sync; %
CALCIUMroiTS.deeplabcut.ligth.x = CSV_file(laserframe,2);
CALCIUMroiTS.deeplabcut.ligth.y = CSV_file(laserframe,3);
CALCIUMroiTS.deeplabcut.tail1.x = CSV_file(laserframe,5);
CALCIUMroiTS.deeplabcut.tail1.y = CSV_file(laserframe,6);
CALCIUMroiTS.deeplabcut.tail2.x = CSV_file(laserframe,8);
CALCIUMroiTS.deeplabcut.tail2.y = CSV_file(laserframe,9);
CALCIUMroiTS.deeplabcut.tail3.x = CSV_file(laserframe,11); %
CALCIUMroiTS.deeplabcut.tail3.y = CSV_file(laserframe,12);
CALCIUMroiTS.deeplabcut.tail4.x = CSV_file(laserframe,14); %
CALCIUMroiTS.deeplabcut.tail4.y = CSV_file(laserframe,15);
CALCIUMroiTS.deeplabcut.tail5.x = CSV_file(laserframe,17); %
CALCIUMroiTS.deeplabcut.tail5.y = CSV_file(laserframe,18);
CALCIUMroiTS.deeplabcut.Eye1.x = CSV_file(laserframe,20); %
CALCIUMroiTS.deeplabcut.Eye1.y = CSV_file(laserframe,21);
CALCIUMroiTS.deeplabcut.Eye2.x = CSV_file(laserframe,23); %
CALCIUMroiTS.deeplabcut.Eye2.y = CSV_file(laserframe,24);
CALCIUMroiTS.deeplabcut.Mouth.x = CSV_file(laserframe,26); %
CALCIUMroiTS.deeplabcut.Mouth.y = CSV_file(laserframe,27);

CALCIUMimg('savewave', CALCIUMroiTS, CALCIUMroiTS.ref,list,nfish);

%3D matrix with coordinates for the ploting

Coor3D = NaN(7,2,length(laserframe)); % 10 points, 2 colums (x,y), length frames
for i = 1:length(laserframe)
    Coor3D(:,1,i) = [CALCIUMroiTS.deeplabcut.Mouth.x(i);CALCIUMroiTS.deeplabcut.ligth.x(i);CALCIUMroiTS.deeplabcut.tail1.x(i);...
        CALCIUMroiTS.deeplabcut.tail2.x(i);CALCIUMroiTS.deeplabcut.tail3.x(i);CALCIUMroiTS.deeplabcut.tail4.x(i) ;...
        CALCIUMroiTS.deeplabcut.tail5.x(i)]; % xcolumn
    Coor3D(:,2,i) = [CALCIUMroiTS.deeplabcut.Mouth.y(i);CALCIUMroiTS.deeplabcut.ligth.y(i);CALCIUMroiTS.deeplabcut.tail1.y(i);...
        CALCIUMroiTS.deeplabcut.tail2.y(i);CALCIUMroiTS.deeplabcut.tail3.y(i);CALCIUMroiTS.deeplabcut.tail4.y(i) ;...
        CALCIUMroiTS.deeplabcut.tail5.y(i)]; % ycolumn
end 

Vec3D = NaN(6,2,length(laserframe)); 
for i = 1:size(Vec3D,1)
   Vec3D(i,1,:) = Coor3D(i,1,:)-Coor3D(i+1,1,:);
   Vec3D(i,2,:) = Coor3D(i,2,:)-Coor3D(i+1,2,:);
end 

 % Using atan2d. ALL CONSECUTIVE BODY PARTS
    Ang3D = NaN(5,1,length(laserframe)); % Has positive and negative angles
    for i = 1:size(Vec3D,1)-1
         Ang3D(i,1,:)= atan2d(Vec3D(i,1,:).*Vec3D(i+1,2,:)-Vec3D(i,2,:).*Vec3D(i+1,1,:),...
             Vec3D(i,1,:).*Vec3D(i+1,1,:)+Vec3D(i,2,:).*Vec3D(i+1,2,:)); 
     end 
    
     % Using atan2d. TWO SELECTED BODY PARTS
     % VECTORS. (i.e. 1--> HEAD-HB  ;  6--> Tail5)
     i = [1 5];
     AngHeadTail3D(1,1,:)= atan2d(Vec3D(i(1),1,:).*Vec3D(i(2),2,:)-Vec3D(i(1),2,:).*Vec3D(i(2),1,:),...
             Vec3D(i(1),1,:).*Vec3D(i(2),1,:)+Vec3D(i(1),2,:).*Vec3D(i(2),2,:)); 
%      plot(squeeze(AngHeadTail3D))
    
     %Calculation of the matrix that contains the SUM of all ANGLE  for each
     %frame
     SumAng3D = sum(Ang3D);
     SumAng3DPosi = sum(abs(Ang3D)); % Has just positive angles
     % size(SumAng3DPosi)
     plot(squeeze(SumAng3D))
     plot(squeeze(AngHeadTail3D))
     plot(1:length(laserframe),CALCIUMroiTS.deeplabcut.tail5.y)
     CALCIUM.SumAng3D = SumAng3D ; 
     CALCIUM. AngHeadTail3D =  AngHeadTail3D; 
     CALCIUM.Ang3D = Ang3D; 

     CALCIUMimg('save',CALCIUM,[],list,nfish);  % If you want to save the raw movie

     %   vidObj = VideoReader(fullfile(Videoreal_directory,a4.name)); % Importing original video
    %    implay(fullfile(Videoreal_directory,a4.name)); 


figure(2)

plot(time_sync,CALCIUMroiTS.deeplabcut.tail3.x ,'r')
title(char(CALCIUM.list(nfish,1)))
ylabel('\Deltax Coordinates tail3 [px]')
xlabel ('Time [s]')
saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),strcat('Tail_Movement_',char(CALCIUM.list(nfish,1)),'.fig')))
close


%% Visualization of the swim trace and the calcium trace 

% Correction problem with the trasposition
if size(CALCIUMroiTS.diff_perc03.data,1)>size(CALCIUMroiTS.diff_perc03.data,2)
    CALCIUMroiTS.diff_perc03.data = CALCIUMroiTS.diff_perc03.data';
    %CALCIUMimg('savewave', CALCIUMroiTS,[],list,nfish);
else
end


%{

Rescaling camera video time.

Because how the experiment is performed:

1. starting video recording
2. starting confocal recording
3. stoping confocal recording
4. stoping video recording

the length of the video should be always longer than 
the lenght of the confocal file and the lenght of the video with the
laser on should be always the same than the length of the confocal file.
Due to an unknow problem (maybe the camera is overheating and can not
keep the rate of adquisition) sometimes this is not the case and the length of
the video with the laser on is shorter than the confocal file. to fi that,
we assumed that the both recordings should be the same length so we make
matching the length of the video with the length of the confocal recording.

This process will create a variable that is called for the camera video
that is called time_scaled. Both now, the camera video withthe laser on time and the
confocal recording should have the same length.

The procedure will be the next one. A varibale CALCIUM.delay will be calculated
next. this varibale gives you the delay btw the confocal recording and
the camera video laser on. For example last time point of the camera video laser on is 10 seconds
and the confocal video is 12 seconds, so the delay will be 2 seconds. this
mean for both of them to be the same that the camera video laser on must be
"streched" an extra of 2 seconds so at the end, both final values are 12
seconds. To do this we divide this delay by the total number of time points
and to each orignal time point we sum this diference in an cumulative way.
For example, if the result of dividng delay/length(orinal_trace) = 0.2 this
means that to the first timepoint of the original trace it will be sum 0.2,
to the second time point it will be added 0.2*2, to the third 0.2*3 and so
on. To the last timepont it will be added 0.2*length(orinal_trace) that is
exatcly the delay that we have measured. This is a way of homogeneously
strecht the trace of the camera video laser on to match the confocal
recording

%}

% This will give us how much delay is the camera laser on video with respect to the confocal recoridng:
% if is <0 the confocal recording is longer than the video that is the
% usual problem
CALCIUM.delay = CALCIUMroiTS.deeplabcut.lasertimesync(end)-CALCIUMroiTS.diff_perc03.times(end);
if CALCIUM.delay<0
else
    error('The camera video laser time on is LONGER than the confocal recording')
end


% Geting the time that must be added cumulativelty
sr_rescaled = -CALCIUM.delay/length((CALCIUMroiTS.deeplabcut.lasertimesync));

% The new variable must have the same amount of time points than the original lasertimesync, so we prelocate with NaN using the same length
time_scaled = NaN(1,length(CALCIUMroiTS.deeplabcut.lasertimesync));
% Geting the first timepoint of the new time variable
time_scaled(1) = sr_rescaled + CALCIUMroiTS.deeplabcut.lasertimesync(1);
% The cumulative variable will be suma
suma =  sr_rescaled;

% Creating the rest of the new time variable
for i = 2:length((CALCIUMroiTS.deeplabcut.lasertimesync))
    suma = suma + sr_rescaled;
    time_scaled(i) = suma + CALCIUMroiTS.deeplabcut.lasertimesync(i) ;
end

% Saving
CALCIUMroiTS.deeplabcut.time_scaled = time_scaled';
CALCIUMimg('savewave', CALCIUMroiTS,[],list,nfish);

% Ploting both the swiming trace (video camera) and the calcium trace (confocal recording)
for i = 1:length(CALCIUMroiTS.roi.labels) % loop acrros the cells

    figure(i)

    % No rescaled time:
    subplot(2,1,2)
    % swim trace
    yyaxis left
    plot(CALCIUMroiTS.deeplabcut.time_scaled,squeeze(CALCIUM.AngHeadTail3D)','g','LineWidth',3); hold on
    yyaxis left
    title(strcat('Calcium activity-Tail movement-SCALED-TIME','-',CALCIUM.ref,'-',CALCIUMroiTS.roi.labels(i)));
    xlabel('Time [s]')
    ylabel('\Deltax Coordinates tail3 [px]')

    yyaxis right % calcium trace
    plot(CALCIUMroiTS.diff_perc03.times,-CALCIUMroiTS.diff_perc03.data(i,:),'b','LineWidth',3);
    yyaxis right
    ylabel('% \Delta F')


    % Rescaled time:
    subplot(2,1,1)
    % swim trace
    yyaxis left
    plot(CALCIUMroiTS.deeplabcut.lasertimesync,squeeze(CALCIUM.AngHeadTail3D)','g','LineWidth',3); hold on
    yyaxis left
    title(strcat('Calcium activity-Tail movement-NON-SCALED-TIME','-',CALCIUM.ref,'-',CALCIUMroiTS.roi.labels(i)));
    xlabel('Time [s]')
    ylabel('\Deltax Coordinates tail3 [px]')

    yyaxis right % calcium trace
    plot(CALCIUMroiTS.diff_perc03.times,-CALCIUMroiTS.diff_perc03.data(i,:),'b','LineWidth',3);
    yyaxis right
    ylabel('% \Delta F')


    saveas(gcf,fullfile(CALCIUM.path.data,'dataCALCIUM',char(CALCIUM.list(nfish,1)),...
        strcat(char(CALCIUMroiTS.roi.labels(i)),'_','Traces',char(CALCIUM.list(nfish,1)),'.fig')))
end

CALCIUMimg('save',CALCIUM,[],list,nfish);


% Export the x variable in CSV file
%  csvwrite(fullfile(CALCIUM.path.data,'dataCALCIUM',...
%      char(CALCIUM.list(nfish,1)),'Tail_x_CSV'),CALCIUMroiTS.deeplabcut.tail3.x);
csvwrite(fullfile(CALCIUM.path.data,'dataCALCIUM',...
    char(CALCIUM.list(nfish,1)),'Tail_CSV'),squeeze(SumAng3D));








