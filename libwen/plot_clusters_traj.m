function plot_clusters_traj(xx,pertab,proj,ci,cm)   
    for k=1:3
        subplot(1,3,k);
        hold on
        uq=unique(ci);
        for i=1:length(uq)
            cli=find(ci==uq(i));
            for j=1:numel(cli)
                idx=pertab(cli(j),1):pertab(cli(j),1)+pertab(cli(j),4);    
                x=xx(idx,:); 
                plot3(x(:,proj(k,1)),x(:,proj(k,2)),x(:,proj(k,3)),'color',cm(i,:));                
            end
        end     
    end
end