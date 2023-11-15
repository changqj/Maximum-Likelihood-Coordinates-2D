function [mvc,thetas] = mvcoordinates(x,p)
% CopyRight:  Qingjun Chang @USI
% meav value coordinates 

n = size(p,2);

% projection step
v = p - repmat(x,1,n);

r0 = vecnorm(v);
v = v./r0;    % projected points on unit circle


thetas = acos(dot(v,v(:,[2:end,1]))).*sign(dot(cross([v;ones(1,n)],...
    [v(:,[2:end,1]);ones(1,n)]),repmat([0;0;1],1,n)));

T = tan(thetas/2);

d = T + T([end 1:end-1]);
d = d./r0;
d = d'/sum(d);
mvc = d;


% d = sin((thetas+thetas([end 1:end-1]))/2);
% d = d./r0;
% d = d'/sum(d);
% mvc = d;

% t(1) = 0.2;
% for i = 2:n
%     t(i) = thetas(i)-t(i-1);
% end
% d = sin(t);
% d = d'/sum(d);

end
