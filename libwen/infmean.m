function mn=infmean(x,dim)
    x(isinf(x))=nan;
    mn=nanmean(x,dim);
end