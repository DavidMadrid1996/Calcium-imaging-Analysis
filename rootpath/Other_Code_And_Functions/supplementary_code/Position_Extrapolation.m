user_settings
data = bfopen(char("C:\Users\Usuario\OneDrive\Escritorio\220825_1. UAS.GCaMP chx10.gal4\Experiment-175.czi"));
% czifinfo(char(strcat(CALCIUM.CZI_directory_loop,'.czi')))

for j = 1:length(data{1,1})
  movies4D (:,:,j) = data{1,1}{j,1};
end 


for j = 1:length(data{1,1})
  movies4D (:,:,j) = data{1,1}{j,1};
end 

movies4D = im2double(movies4D);
Outputzstack = NaN(size(movies4D,1),size(movies4D,2),size(movies4D,3)); 
for i = 1:size(movies4D,3)
Outputzstack(:,:,i) = imadjust(movies4D(:,:,i),[0.0000 0.1847 ],[]);
end 
% imshow(Outputzstack)
% imcontrast

clear

snap = bfopen(char("C:\Users\Usuario\OneDrive\Escritorio\220825_1. UAS.GCaMP chx10.gal4\Snap-26.czi"));
% czifinfo(char(strcat(CALCIUM.CZI_directory_loop,'.czi')))

snap4D = snap{1,1}{1,1};

Outputsnap= imadjust(snap4D,[0.0000 0.1847 ],[]);
imshow(Outputsnap)



