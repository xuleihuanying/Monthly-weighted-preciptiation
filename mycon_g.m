% Here is the constraint function as Eq. (9)
function [c,ceq] = mycon_g(r,s)
% s, covariance, (N-1)*(N-1)
% r, 1*N, free parameters
m = length(r);
u = ones(m-1, 1);
k = (det(s))^(1 / (m-1) );
c = -( r(end) - ((r(1:end-1))' - r(end)*u)' * inv(s) * ( (r(1:end-1))' - r(end)*u) ) / k;
ceq = [];
end