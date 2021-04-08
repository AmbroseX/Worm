close all
clear segmented_curv t_w_segmented_curv t_w_DLPon_frame
segment_length=20;
ant_boundary = 40;
pos_boundary = 60;
start_segment=1:segment_length:100;

Num_segments=length(start_segment);


for k = 1:cyclenum
    if t_w_illu_time(k,3) == ant_boundary && t_w_illu_time(k,4) == pos_boundary
        segmented_curv=zeros(size(t_w_curvdatafiltered{k},1),Num_segments);
        for i=1:Num_segments
            if ~isempty(t_w_curvdatafiltered{k})
            for j=1:size(t_w_curvdatafiltered{k},1)

                    segmented_curv(j,i)=mean(t_w_curvdatafiltered{k}(j,start_segment(i):start_segment(i)+segment_length-1));

            end
            end
        end
    DLPon_frame = find(time_auto{k} == 0);



    h1=figure;

    plot(1:j,segmented_curv(:,1),'k');
    hold on
    plot(1:j,segmented_curv(:,2),'r');
    plot(1:j,segmented_curv(:,3),'b');
    plot(1:j,segmented_curv(:,4),'g');
    plot(1:j,segmented_curv(:,5),'c');
    legend('0-20','20-40','40-60','60-80','80-100');
    plot([DLPon_frame,DLPon_frame],[-8,10],'r','Linewidth',2);
    t_w_segmented_curv{k} = segmented_curv;
    t_w_DLPon_frame(k) = DLPon_frame;
    savefig(h1,strcat(wormname,'_',num2str(k),' seg_cur'));
    saveas(h1,strcat(wormname,'_',num2str(k),' seg_cur'),'jpg');
    end
end
save([wormname,'seg_curv'],'t_w_segmented_curv','t_w_DLPon_frame');
