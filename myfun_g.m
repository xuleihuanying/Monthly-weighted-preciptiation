% Here is the objective function as Eq. (8)
function F = myfun_g(r,s)
% s, covariance, (N-1)*(N-1)
% r, 1*N, free parameters
m = length(r);
u = ones(m-1, 1);
r_2 = s - r(end)*(u*u') + u*(r(1:end-1)) + (r(1:end-1))'*u';
r_all = [r_2 (r(1:end-1))'; r(1:end-1) r(end) ];
ff = 0;
k = (det(s))^(1 / (m-1) );
for i=1:m-1
    for j=i+1:m
        ff = ff + r_all(i,j) * r_all(i,j) / k/k;
    end
end
F = ff;

end