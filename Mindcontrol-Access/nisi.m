%close all
filename_simplicity=wormname;
save(filename_simplicity);

evalc([wormname,'_amp=t_w_amp']);
evalc([wormname,'_curdata=t_w_curdata']);
evalc([wormname,'_curvdatafiltered=t_w_curvdatafiltered']);
%evalc([wormname,'_angle_output=t_w_angledata']);
if exist('t_w_illu_time','var')
    evalc([wormname,'_illu=t_w_illu_time']);
    
    save([wormname,'_amp&cur'],[wormname,'_amp'],[wormname,'_curdata'],[wormname,'_curvdatafiltered'],[wormname,'_illu']);
elseif exist('iteration_info','var')
    evalc([wormname,'_iteration=iteration_info']);
    save([wormname,'_amp&cur'],[wormname,'_amp'],[wormname,'_curdata'],[wormname,'_curvdatafiltered'],[wormname,'_iteration_info']);

    
end
clear ans



