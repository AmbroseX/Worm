function y=asort(y,x,asc)
    %sort the matrix y according to columns of x
    if ~exist('asc','var')
        asc='ascend';
    else
        asc='descend';
    end
    [~,si]=sort(x,asc);
    y=y(si,:);
end