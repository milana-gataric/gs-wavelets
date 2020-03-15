%Calculates pointwise evluations of the 2D Fourier transform at the points
%which can be nonuniform

%%NB: this works for most simple functions and for a smaller number of sampling
%     points, but may sometimes return NaN values,
%     So use this at your own peril... 

%Inputs:
%   sp:     2D vector of sampling points (where the FT is calculated)
%   f_expr: this is a string and should contain only 2 variables, x, y. e.g. 'x+1', 'x*y'
%   a1, b1: interval on which Fourier transform is evaluated on in x
%   a2, b2: interval on which Fourier transform is evaluated on in y
%


function omega = Samples_NU(sp, f_expr, a1, b1, a2, b2 )

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

        omega=FT(2*pi*sp(:,1), 2*pi*sp(:,2));
        w=2*pi*sp;
        omega(sp(:,2)==0)=FT_y0(w(sp(:,2)==0,1));
        omega(sp(:,1)==0) = FT_x0(w(sp(:,1)==0,2));
        omega(sp(:,1)==0&sp(:,2)==0)=FT_expr_00;        


end

% Copyright (c) 2014. Clarice Poon and Milana Gataric
