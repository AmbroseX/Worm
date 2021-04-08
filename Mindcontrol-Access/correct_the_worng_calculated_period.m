[se,status,pklcs]=draw_number_peaks;
answer=inputdlg({'do you konw the start valid start peakpoint number?','and the end peakpoint number','any false point?(optional)'}...
    ,'',1);
pk_start=str2num(answer{1});
pk_end=str2num(answer{2});
false_point=str2num(answer{3});

fprintf('the %d th corrected frequency is: ',se);
if status==1
    t_w_meanperiod_DLPoff(se,1)=(pklcs(pk_end,2)-pklcs(pk_start,2))/(pk_end-pk_start-length(false_point))/fps;
    t_w_meanperiod_DLPoff(se,2)=1/t_w_meanperiod_DLPoff(se,1);
    display(t_w_meanperiod_DLPoff(se,2));
elseif status==3 || status==2
    t_w_meanperiod_DLPon(se,1)=(pklcs(pk_end,2)-pklcs(pk_start,2))/((pk_end-pk_start-length(false_point))*fps);   
    t_w_meanperiod_DLPon(se,2)=1/t_w_meanperiod_DLPon(se,1);
    display(t_w_meanperiod_DLPon(se,2));
    t_w_mat(se,6) = 1/t_w_meanperiod_DLPon(se,1);
end

close