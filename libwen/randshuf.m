function shuf=randshuf(x)
    sz = size(x);
    for i=1:sz(2)
        randi = randperm(sz(1));
        shuf(:,i) = x(randi,i);
    end
    %randi = randperm(sz(1)*sz(2));
    %shuf = x(randi);
    %shuf = reshape(shuf,sz(1),sz(2));
%     plot(shuf)
%     histogram(x(:,1))
end