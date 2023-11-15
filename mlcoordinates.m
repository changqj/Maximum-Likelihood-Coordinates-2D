function [mlc] = mlcoordinates(x,v,distances)
% CopyRight:  Qingjun Chang @USI
% maximum likelihood coordinates 

n = size(v,2);
d = size(v,1);
idx = @(i) mod(i-1,n)+1;
v_cir = v - repmat(x,1,n);
r = vecnorm(v_cir);     % the distances from x to each v_i
v_cir = v_cir./r;     % the vertices of the projected polygon

l = v_cir(:,[n 1:n-1]) + v_cir;
vecnorml = vecnorm(l);
s = l./vecnorml;
t = s(:,[2:n 1]);

v_bar = s+t;
length_s_plus_t = vecnorm(v_bar);
v_bar = v_bar./length_s_plus_t;     % the vertices of the polygon after smooth step


v_hat = v_bar./distances;


F = @(phi) -sum(log(n+[phi(1) phi(2)]*v_hat));

% You can also use the following Newton's method code instead
options = optimoptions("fminunc",OptimalityTolerance=1e-10,Display="off");
[phi,~,~,~] = fminunc(F,[0,0],options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k = 0;
% 
% error = 1e-8;
% phi = zeros(1,d);
% 
% F1 = @(phi) -(sum(v_hat./(n+phi*v_hat),2));   % d-by-1
% F2 = @(lambda,i,j) sum(v_hat(i,:).*v_hat(j,:)./((n+phi*v_hat).^2));    % d-by-d
% 
% H = zeros(d,d);
% while 1
% 
%     g = F1(phi);
%     k = k + 1;
%     E(k) = norm(g);
%     if E(k)<=error
%         break;
%     end
%     for i = 1:d
%         for j = 1:d
%             H(i,j) = F2(phi,i,j);
%         end
%     end
%     Delta = -H\g;
% 
%     % Armijo linear search
%     mm = 0; mk = 0;
%     rho = 0.55; sigma = 0.4;
%     while(mm<20)
%         if(F(phi+rho^mm*Delta')<F(phi)+sigma*rho^mm*g'*Delta)
%             mk = mm;
%             break;
%         end
%         mm = mm +1;
%     end
%     alpha = rho^mk;
%     phi = phi + alpha*Delta';
% 
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi = phi';
lambda_hat = 1./(n+phi'*v_hat);

lambda_bar = lambda_hat./distances;

X = 1./(vecnorml.*length_s_plus_t([n 1:n-1]));
Y = (1./vecnorml + 1./vecnorml([2:n 1]))  ./ length_s_plus_t;
Z = 1./(vecnorml([2:n 1]).*length_s_plus_t([2:n 1]));
lambda_cir = lambda_bar;
for i = 1:n
    lambda_cir(i) = lambda_bar(idx(i-1))*X(i)+lambda_bar(i)*Y(i)+lambda_bar(idx(i+1))*Z(i);
end
lambda = lambda_cir./r;

mlc = lambda/sum(lambda);

end



