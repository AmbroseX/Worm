if ~exist('mcd','var')
    wormspeedanswer = inputdlg({'which video do you want to load'}, '', 1);
    
    wormname = [wormspeedanswer{1},'.mat'];
    load(wormname);
    clearvars -except mcd t_w_illu_time worm_length wormspeedanswer fname;
    wormcheckpoints(1,:) = t_w_illu_time(:,1)-120;
    wormcheckpoints(2,:) = t_w_illu_time(:,1)+200;
    wormcheckpoints(3,:) = t_w_illu_time(:,1)+300;
    wormcheckpoints(5,:) = t_w_illu_time(:,1)-350;
    wormcheckpoints(6,:) = t_w_illu_time(:,1)+350;
    
end