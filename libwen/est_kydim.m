function [dky_lin,dky_cub,md]=est_kydim(lap)
    rng=1:1e-4:size(lap,2);
    dky_lin=interp1(1:size(lap,2),cumsum(lap),rng,'linear');
    dky_cub=interp1(1:size(lap,2),cumsum(lap),rng,'pchip');
    dky_lin=rng(find(dky_lin<0,1,'first'));
    dky_cub=rng(find(dky_cub<0,1,'first'));
    if isempty(dky_lin),dky_lin=nan;end
    if isempty(dky_cub),dky_cub=nan;end
    md=amax(cumsum(lap));
end
