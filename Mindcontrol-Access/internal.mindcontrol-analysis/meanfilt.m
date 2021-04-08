%Go through a set of data and take the mean for a given window size

%After you export the analysis from the previewr program
%load the .mat file
%and run this:
q=handles.curvdata;
w=25;%window size
N=size(q,1); %number of frames
for k=1:N-w;
    a(k)=mean(mean(q(k:k+w,:)));
end
plot(a);