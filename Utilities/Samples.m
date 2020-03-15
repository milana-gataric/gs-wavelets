%Calculates pointwise evluations of the Fourier transform
%Inputs:
%   eps:    Fourier sampling density
%   M:      Range of Fourier samples
%   f_expr: this is a string and should contain only 2 variables, x, y. e.g. 'x+1', 'x*y'
%   a1, b1: interval on which Fourier transform is evaluated on in x
%   a2, b2: interval on which Fourier transform is evaluated on in y
%
%Examples: 
% Samples( eps, M, 'x+1', a, b ) returns a vector, for k=M:1:-M+1, the
% entries are:
%      sqrt(eps) *  int_a^b (x+1) exp(-2*pi*eps*k*x) dx .
%
% Samples( eps, M, '(1-x)*(1-y)', a1, b1, a2, b2 ) returns a matrix, 
% for k1, k2 =M:1:-M+1, the (k1,k2) entry is
%   eps * int_a2^b2 int_a1^b1 (1-x)*(1-y) exp(-2*pi*eps*(k1*x+k2*y)) dx
%
function omega = Samples( eps, M, f_expr, a1, b1, a2, b2 )

    if nargin < 7

        f=str2func( strcat('@(x)', f_expr) );
        syms x z;
        FTexpr = char(int(f(x)* exp(-1i*x*z), x, a1, b1));
        FT_0 = int(f(x), x, a1, b1);

        FTexpr = strrep(FTexpr, '*', '.*');
        FTexpr = strrep(FTexpr, '^', '.^');
        FTexpr = strrep(FTexpr, '/', './');



        FT=str2func( strcat('@(z)', FTexpr) );
        omega = sqrt(eps)*FT(2*pi*eps*(M:-1:-M+1)');
        omega(M+1) = sqrt(eps)*FT_0;
    else

        f=str2func(strcat('@(x,y)', f_expr));
        syms x y z1 z2;
        FT_expr = char(int( int(f(x,y)*exp(-1i*x*z1)*exp(-1i*y*z2), x, a1, b1), y, a2, b2));
        FT_expr_x0 = char(int( int(f(x,y)*exp(-1i*y*z2), x, a1, b1), y, a2, b2));
        FT_expr_y0 = char(int( int(f(x,y)*exp(-1i*x*z1), x, a1, b1), y, a2, b2));
        FT_expr_00 = int( int(f(x,y), x, a1, b1), y, a2, b2);


        FT_expr = strrep(FT_expr, '*', '.*');
        FT_expr = strrep(FT_expr, '^', '.^');
        FT_expr = strrep(FT_expr, '/', './');
        FT_expr_x0 = strrep(FT_expr_x0, '*', '.*');
        FT_expr_x0 = strrep(FT_expr_x0, '^', '.^');
        FT_expr_x0 = strrep(FT_expr_x0, '/', './');
        FT_expr_y0 = strrep(FT_expr_y0, '*', '.*');
        FT_expr_y0 = strrep(FT_expr_y0, '^', '.^');
        FT_expr_y0 = strrep(FT_expr_y0, '/', './');

        FT = str2func(strcat( '@(z1, z2)', FT_expr));
        FT_y0 = str2func(strcat( '@(z1)', FT_expr_y0));
        FT_x0 = str2func(strcat( '@(z2)', FT_expr_x0));


        omega = zeros(2*M, 2*M);
        for k=M:-1:-M+1
            omega(:,M-k+1) = FT(2*pi*eps*(M:-1:-M+1)', 2*k*eps*pi);
        end
        w=2*pi*eps*(M:-1:-M+1)';

        omega(:,M+1)=FT_y0(w);
        omega(M+1,:) = FT_x0(w);
        omega(M+1,M+1)=FT_expr_00;
        omega = eps*omega;

    end
end


% Copyright (c) 2014. Clarice Poon and Milana Gataric
