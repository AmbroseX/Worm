function segdata=savesegdata(j,j1,segmented_curv,time,sequence)
segdata=cell(20,1);
segdata_1=zeros(j,5);

for j_length=1:j
    
if time(j_length)>-3.02
for portion_worm=1:5
    for k=1:(j-j1)
        %fprintf('%d\n',segmented_curv(j_length,portion_worm));
    segdata_1(k,portion_worm)=segmented_curv(j_length,portion_worm);
   end
    
segdata{sequence,1}=segdata_1;
end
end
end