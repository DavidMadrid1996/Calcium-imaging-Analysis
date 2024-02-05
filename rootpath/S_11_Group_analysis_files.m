%% Data Extraction Script. Last modified 2024/01/17

%{

                                    INFO:  

 This script will load the names and path to the files that you want to analyse.
 For the script to run it is required having the list of animals that you
 want to analysis.There will anatomatically load the list of animals in different conditions 
 (sw, turning, etc. even if it will use just one of the conditions), and there will find the 
 path to each of this animals, trial and swimming episode.


%}

%% 
clear
user_settings
% Here you must set the list of animals that you want to analyzed
load('Possible_candidates_sw.mat');
load('Possible_candidates_struggle.mat');
load('Possible_candidates_turning_left.mat')
load('Possible_candidates_turning_right.mat')
load('temp_list.mat')
load('temp_list_ex_vivo.mat')
load('Just_one_episode_Possible_candidates_sw.mat')

% Here you select just the one that it will be look for inside all the
% folders. The File_Location(XXX, '\**\*.mat'); XXX must be one of the
% animals list:
[ValidPathdataCalcium, ValidPathdatawaves,Animal_List_unique] = File_Location(temp_list_ex_vivo, '\**\*.mat');

% Output of this code will be:

% * ValidPathdataCalcium --> path to the Calcium matlab file for each
%                            animal -trial-swimming epsiode

% * ValidPathdatawaves -->   path to the Wave matlab file for each 
%                            animal -trial-swimming epsiode

% * Animal_List_unique --> New list that does not include repetitions. This
%                          is, it will include just the 

