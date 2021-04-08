if ~exist('mcd','var')
    answer = inputdlg({'which video do you want to load'}, '', 1);
    
    wormname = [answer{1},'.mat'];
    load(wormname);
    clearvars -except mcd t_w_illu_time worm_length
    istarts = input('istarts\n');
    iends = istarts+400;
end
    spd = cell(1,length(istarts));
for counting = 1:length(istarts)
    istart = istarts(counting);
    iend = iends(counting);
    stagepos = worm_stagepos{counting};
    TrackWormRealSpace_true;
    spd{counting} = speed;
    figure;plot(speed);
end