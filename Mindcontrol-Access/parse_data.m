function [g1,g2]=parse_data(group1,group2)

k1=length(find(~isnan(group1)));
k2=length(find(~isnan(group2)));
k=min(k1,k2);
g1=zeros(k,1);
g2=g1;
j=0;
for i=1:length(group1)
    if (~isnan(group1(i)))&&(~isnan(group2(i)))
        j=j+1;
        g1(j)=group1(i);
        g2(j)=group2(i);
    end
end


    
