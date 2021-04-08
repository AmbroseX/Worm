function [Ratio, RatioTimeStamps]=TimeSeriesBinaryEvents2Ratio(T,R,nbins)
% [Ratio, RatioTimeStamps]=TimeSeriesBinaryEvents2Ratio(t,r,nbins)
% 
% Convert a time series of binary events into a binned series of ratios of 
% events to non events.
%
% r is a vector of 0 or 1 binary events
% t is a vector contaning time corresponding to the 0 or 1 binary events
% nbins is the number of bins desired
%
% Andrew Leifer
% leifer@fas.harvard.edu
nb=nbins; %number of time bins

%Sort the responses as a function of time
[Tsort Ix]=sort(T);
Rsort=R(Ix);

RsortSum=cumsum(Rsort);
m=0;
for k=1:nb
    %These indices correspond to the location of every 200seconds in Tsort
    TsortIxEven(k)=findClosest(Tsort,k*Tsort(end)/nb);
end
NumEventsPerBin=RsortSum(TsortIxEven)-RsortSum([1 TsortIxEven(1:end-1)]);
NumStimuliPerBin=TsortIxEven - [1 TsortIxEven(1:end-1)] ;
Ratio=NumEventsPerBin./(NumStimuliPerBin);
RatioTimeStamps=  ( Tsort(TsortIxEven)+Tsort([1 TsortIxEven(1:end-1)]) )./2;