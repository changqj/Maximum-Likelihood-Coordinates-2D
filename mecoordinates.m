function [w] = mecoordinates(x,p,priorfun_option)
% CopyRight:  Qingjun Chang @USI
% maximum entropy coordinates with prior function
%   x: the point needed to compute coordiantes
%   p: points
%   priorfun_option: the index of prior function to use
%       #1: prior function MEC-1
%       #2: prior function MEC-2
%       #3: prior function MEC-constant
%       #4: prior function MEC-Gaussian


n = size(p,2);
d = size(p,1);

switch priorfun_option
    case 1
        % prior function  MEC-1
        priorfun = vecnorm(repmat(x,1,n)-p)+vecnorm(repmat(x,1,n)-p(:,[2:n 1]))-vecnorm(p(:,[2:n 1]) - p);
        pi_ = 1 ./ (priorfun([n 1:n-1]) .* priorfun);
        m = pi_ / sum(pi_);    % 1-by-n
    case 2
        % prior function  MEC-2
        priorfun = vecnorm(repmat(x,1,n)-p) .* vecnorm(repmat(x,1,n)-p(:,[2:n 1])) + ...
            dot((repmat(x,1,n)-p) ,  (repmat(x,1,n) - p(:,[2:n 1])));
        pi_ = 1 ./ (priorfun([n 1:n-1]) .* priorfun);
        m = pi_ / sum(pi_);    % 1-by-n
    case 3
        % prior function  MEC-constant
        m = ones(1,n);
    case 4
        % prior function  MEC-Gaussian
        beta = 1;
        m = exp(-beta*(vecnorm(p-repmat(x,1,n)).^2));
    otherwise
        m = ones(1,n);
end




v = p - repmat(x,1,n);

k = 0;

error = 1e-8;
lambda = zeros(1,d);
Zi = @(lambda) m .* exp(-lambda*v);     % 1-by-n
Z = @(lambda) sum(Zi(lambda));     % number
F = @(lambda) log(Z(lambda));
F1 = @(lambda) -(v*(Zi(lambda))')/Z(lambda);   % d-by-1
F2 = @(lambda,i,j) (sum(v(i,:).*v(j,:).*Zi(lambda))*Z(lambda)-...
    (v(i,:)*(Zi(lambda))')*(v(j,:)*(Zi(lambda))'))/((Z(lambda))^2);    % d-by-d

H = zeros(d,d);
while 1

    g = F1(lambda);
    k = k + 1;
    E(k) = norm(g);
    if E(k)<=error
        break;
    end
    for i = 1:d
        for j = 1:d
            H(i,j) = F2(lambda,i,j);
        end
    end
    Delta = -H\g;

    % Armijo linear search
    mm = 0; mk = 0;
    rho = 0.55; sigma = 0.4;
    while(mm<20)
        if(F(lambda+rho^mm*Delta')<F(lambda)+sigma*rho^mm*g'*Delta)
            mk = mm;
            break;
        end
        mm = mm +1;
    end
    alpha = rho^mk;
    lambda = lambda + alpha*Delta';

end

w = Zi(lambda)'/Z(lambda);


w = w/sum(w);
end