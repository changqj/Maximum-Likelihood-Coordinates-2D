function [ic] = icoordinates(x,v,K)
% CopyRight:  Qingjun Chang @USI
% iterative coordinates 



n = size(v,2);
v_cir = v - repmat(x,1,n);
r = vecnorm(v_cir);     % the distances from x to each v_i
v_cir = v_cir./r;     % the vertices of the projected polygon

r1 = zeros(K,n);
for i = 1:K
    v_cir = v_cir(:,[2:n 1]) + v_cir;
    r1(i,:) = vecnorm(v_cir);
    v_cir = v_cir./r1(i,:);
end

w = mvcoordinates([0;0],v_cir);

for i = K:-1:1
    w = w./r1(i,:)'/2;
    w = w + w([end 1:end-1]);
end

w = w./r';

ic = w/sum(w);

end



