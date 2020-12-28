% Largest EigenValue CDF of Wishart Matrix
function cdf = CDF_max(x)

global m n;


funcky = @(a,b)  gammainc(x,n-m+a+b-1,'lower') .* gamma(n-m+a+b-1);

divider = 1;

for kk = 1:1:m

    divider = divider * gamma(n-kk+1)* gamma(m-kk+1);

end

numerator_matrix = zeros(m,m);

for ii = 1:1:m

    for jj = 1:1:m


        numerator_matrix(ii,jj) = funcky(ii,jj);


    end

end

cdf = det(numerator_matrix) ./ divider; % From Tall's