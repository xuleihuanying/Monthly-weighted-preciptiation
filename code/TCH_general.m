%
% This code is the implementation of generalized three-cornered hat (TCH) method
% to obtain the random error variance of multiple precipitation products or
% other variables

% input: 
%       x, [m n], n kinds of data with m length
% output:
%       R, [n n], variance-covariance matrix of the error of x
function [S,R] = TCH_general(x) 
% The "reference" can be anyone of the four inputs, here we choose the
% last one
[m1 m2] = size(x);
y = NaN(m1, m2-1);
for i=1:m2-1
    y(:, i) = x(:,i) - x(:,end);
end

S = cov(y);
u = ones(m2-1, 1);

x0 = zeros(1, m2);
x0(end) = (2*u'*inv(S)*u)^-1;

% Here we use the matlab optmization toolbox for computing the N-free
% parameters.
opts = optimset('Algorithm','active-set','TolX',1e-10,'TolCon',1e-10,'Display','off');
x2 = fmincon(@(r) myfun_g(r,S),x0,[],[],[],[],[],[],@(r) mycon_g(r,S),opts);
% Once the free-parameters have been determined, we can compute the remain
% elements of R using
r = x2;
r_2 = S - r(end)*(u*u') + u*(r(1:end-1)) + (r(1:end-1))'*u';
r_all = [r_2 (r(1:end-1))'; r(1:end-1) r(end) ];
R = r_all;

end