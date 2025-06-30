function v = chi_mandelbrot(x, y, M)
    n = size(M,1);
    xi = x * (n-1) + 1;
    yi = y * (n-1) + 1;
    i  = floor(xi);
    j  = floor(yi);
    i(i < 1) = 1;
    i(i > n - 1) = n - 1;
    j(j < 1) = 1;
    j(j > n - 1) = n - 1;

    dx = xi - i;
    dy = yi - j;

    v00 = M( sub2ind([n n], i    , j    ) );
    v10 = M( sub2ind([n n], i + 1, j    ) );
    v01 = M( sub2ind([n n], i    , j + 1) );
    v11 = M( sub2ind([n n], i + 1, j + 1) );

    v = (1-dx).*(1-dy).*v00 + dx.*(1-dy).*v10 + (1-dx).*dy.*v01 + dx.*dy.*v11;
end
