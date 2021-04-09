function x=norm2one(x)
%     row=size(x,1);
%     col=size(x,2);
%     for i=1:row
%         ma=max(x(i,:));
%         mi=min(x(i,:));
%         x(i,:)=(x(i,:) - mi) / (ma - mi);
%     end
ma=max(x);
mi=min(x);

x=(x-mi)./(ma-mi);
end