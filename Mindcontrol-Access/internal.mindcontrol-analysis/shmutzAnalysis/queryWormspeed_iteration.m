LoadWormdata;
close all
if ~exist('wromspeed','var')
    wormspeed = cell(1);
    worm_stagepos = cell(1);
    running_count = 1;
end
answers = inputdlg({'what are numbers of trails do you want to plot?',...    
    'method use: stagepos in yaml (1) or reference point (2)',...
    'Calculation Method: using multiple body points (segmented) (1) or centriod (2)',...
    'save? 1 for yes otherwise no'},'',1);
trailnums = str2num(answers{1});
reference_point_mtd = str2num(answers{2});
CalculationMethodSelection = str2num(answers{3});
for trailnum = 1:length(trailnums)
    switch reference_point_mtd
        case 1
            [speed,stagepos] = queryWormspeed(mcd,worm_length,reference_point_mtd,CalculationMethodSelection,...
            wormcheckpoints(5,trailnums(trailnum)),wormcheckpoints(6,trailnums(trailnum)));
        case 2
            [speed,stagepos] = queryWormspeed(mcd,worm_length,reference_point_mtd,CalculationMethodSelection,...
            wormcheckpoints(1,trailnums(trailnum)),wormcheckpoints(2,trailnums(trailnum)),wormcheckpoints(3,trailnums(trailnum)),fname);
    end
    wormspeed{running_count} = speed;
    worm_stagepos{running_count} =stagepos;
    running_count = running_count + 1;
end
evalc([wormspeedanswer{1},'wormspeed = wormspeed']);            
disp('# of trail polted');
disp(trailnums);
disp('corresponding sequence');
disp([1:trailnum]);
if answers{4} == '1'
    save([wormspeedanswer{1},'speed'],'wormspeed',[wormspeedanswer{1},'wormspeed']);
end