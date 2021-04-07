workpath=fullfile('G:','Data','WenLab','Worm_Embed');
addpath(genpath(fullfile(workpath,'libwen')));
filepath='testDLP';
pathname=fullfile(workpath,'rawdata',filepath); %the rawdata's path
yamlfiles = dir(fullfile(pathname,'*.yaml'));


start_yaml = input('From which file start # of *.yaml file:'); %you can input the star file number
if isempty(start_yaml)  %if the input is empty then start from the first file
    start_yaml = 1;
end