function [demb,orig,orig_idx,sten]=embed_cells(data,wind)   
%[demb,orig,orig_idx,sten]=embed_cells(data,wind) 
%
% Wrapper on delayembed.m for multiexperiment data (with NaNs as missing
% data). This is main function for delay embedding (tau is set to 1 by
% default). 
% 
% INPUT:
%   data:   Ne x 1 cell array of Ne experiments. If there's only 1 experiment then data can be a matrix, the code will convert it to 1x1 cell array.
%           Each of the Ne matrices are of the form  T_i x d, for T_i time
%           points and d measurements. 
%   wind: The embedding window K.
%        
% OUTPUT:
%   demb:     Delay matrix concatenated across all experiments
%   orig:     The original measurements corresponding to the delays (the first
%             d columns of the delay matrix.
%   orig_idx: The corresponding indices in the observation matrix
%   sten:     Starting and ending indices of each experiment.
%
    if ~iscell(data)
        data={data};
    end   
    nd=size(data{1},2)+1;
    [B,imat]=nanSegment(data,nd,wind); %break at each nan segment and concatenate them across all items 
    [demb,orig,sten]=embd(B,wind,nd,length(data),imat); %delay embed observations in B and concatenate into one big matrix
    orig_idx=demb(:,1);
    demb(:,1:(nd):size(demb,2))=[];
end

function [demb,orig,sten]=embd(B,wind,nd,dlen,imat)

    demb=[]; %contains delay embedded data    
    orig=[]; %contains true copy (shifted by delay length)    
    
    len=zeros(length(B),1);
    for scnt=1:length(B)
        seg=B{scnt};
        %perform delay embedding
        Z = delayembed(seg, 1, wind);        
        demb = [demb;Z];
        orig = [orig;Z(:,2:nd)];
        len(scnt)=size(Z,1);
    end
    
    for wcnt=1:dlen
        newlen(wcnt)=sum(len(imat==wcnt));
    end
  
    sten=[circshift(cumsum(newlen'),1)+1 cumsum(newlen')];
    sten(1)=1;

end

function [B,imat]=nanSegment(data,nd,wind)
    B={}; 
    imat=[];    
    
    
    for wcnt=1:length(data)
        th=data{wcnt}(:,1:nd-1);
        naninds=find(isnan(th(:,1)))'; %find nan indices
        th=[(1:size(th,1))' th];
        th(naninds,1)=nan;
        if ~isempty(naninds)
            if naninds(1)~=1 %if the first element is not nan
                nanbreaks=find(diff([1 naninds]')~=1); %find nan break
            else
                nanbreaks=[1;find(diff([naninds]')~=1)+1]; %find nan break and include the 1st index
            end
            %collect starting ending indices of sequence of nan's
            sten=[naninds(nanbreaks)' naninds(circshift(nanbreaks-1,-1)+[zeros(length(nanbreaks)-1,1);1])'];
            sten(end)=naninds(length(naninds));
            
            %make a new celll array for each nan segment
            for i=1:size(sten,1)+1
                if i==1                   
                    A{i}=th(1:sten(i,1)-1,1:nd); 
                elseif i==size(sten,1)+1
                    if sten(i-1,2)<length(th)                      
                        A{i}=th(sten(i-1,2)+1:end,1:nd);
                    end
                else                    
                    A{i}=th(sten(i-1,2)+1:sten(i,1)-1,1:nd);
                end
            end
            if isempty(A{1}),A(1)=[];end
        else
            A{1}=th;
        end
       
        B=[B;A'];
        imat=[imat;repmat(wcnt,length(A),1)];
        imat=imat(cellfun(@(x) size(x,1),B)>25);
        B=B(cellfun(@(x) size(x,1),B)>25);
        clear A;
    end
end