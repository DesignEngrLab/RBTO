function [roots,weights] = GausQuadBasics(p,type,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% returns the roots and weights needed to implement Gaussian quadrature 
%% corresponding to the orthogonal polynomials 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% p - degree of the polynomial used for the approximation
%%
%% type - type of polynomial
%% 'Hermite': weight function is w(x)=exp(-x^2/2)/(2pi)^.5
%% 'Legendre': weight function is w(x)=1
%% 
%% in the case of type = 'Legendre' it's possible to change the interval of
%% integration to [a,b]
%%
%% varargin - optional list of arguments 
%% contains the endpoints [a,b] of the interval of integration 
%% for type = 'Legendre'
%% if varargin is empty and type = 'Legendre', [a,b] is assumed [-1,1]
%%
%% examples of run:
%% [roots,weights] = GausQuadBasics(3,'legendre');
%% [roots,weights] = GausQuadBasics(2,'hermite');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin == 2)&& strcmpi(type,'Legendre')
    a = -1; b = 1;
elseif (nargin > 2)
    a = varargin{1};
    b = varargin{2};
end

roots = zeros(p,1);
weights = zeros(p,1);
switch lower(type)
    case 'legendre'
        switch p
            case 1
                roots(1) = 0; weights(1) = 2;
            case 2
                roots(1) = 1/sqrt(3); weights(1) = 1;
                roots(2) = -roots(1); weights(2) = weights(1);
            case 3                
                roots(1) = sqrt(3/5); weights(1) = 5/9;
                roots(2) = -roots(1); weights(2) = weights(1);
                roots(3) = 0;         weights(3) = 8/9;
            case 4
                roots(1) = sqrt((3-2*sqrt(6/5))/7); weights(1) = (18+sqrt(30))/36;
                roots(2) = -roots(1);               weights(2) = weights(1);
                roots(3) = sqrt((3+2*sqrt(6/5))/7); weights(3) = (18-sqrt(30))/36;
                roots(4) = -roots(3);               weights(4) = weights(3);
            case 5
                roots(1) = sqrt(5-2*sqrt(10/7))/3;  weights(1) = (322+13*sqrt(70))/900;
                roots(2) = -roots(1);               weights(2) = weights(1);
                roots(3) = sqrt(5+2*sqrt(10/7))/3;  weights(3) = (322-13*sqrt(70))/900;
                roots(4) = -roots(3);               weights(4) = weights(3);
                roots(5) = 0;                       weights(5) = 128/225;
            otherwise
                u = 1:p-1;
                u = u./sqrt(4*u.^2 - 1);

                A = zeros(p,p);
                A(2:p+1:p*(p-1)) = u;
                A(p+1:p+1:p^2-1) = u;

                [v,roots] = eig(A);
                [roots,k] = sort(diag(roots));
                weights = 2*v(1,k)'.^2;
                
        end
         roots = 0.5*((b-a)*roots+(b+a));
         weights = 0.5*weights;
        
    case 'hermite'
        switch p
            case 1
                roots(1) = 0; weights(1) = 1;
            case 2
                roots(1) = 1;         weights(1) = 0.5;
                roots(2) = -roots(1); weights(2) = weights(1);
            case 3                
                roots(1) = sqrt(3);   weights(1) = 1/6;
                roots(2) = -roots(1); weights(2) = weights(1);
                roots(3) = 0;         weights(3) = 2/3;
            case 4
                roots(1) = sqrt((3-sqrt(6))); weights(1) = 1/(4*(3-sqrt(6)));
                roots(2) = -roots(1);           weights(2) = weights(1);
                roots(3) = sqrt((3+sqrt(6))); weights(3) = 1/(4*(3+sqrt(6)));
                roots(4) = -roots(3);           weights(4) = weights(3);        
            case 5
                roots(1) = 0.958572464613819*sqrt(2);  weights(1) = 0.3936193231522/sqrt(pi);
                roots(2) = -roots(1);          weights(2) = weights(1);
                roots(3) = 2.020182870456086*sqrt(2);  weights(3) = 0.01995324205905/sqrt(pi);
                roots(4) = -roots(3);          weights(4) = weights(3);
                roots(5) = 0;                  weights(5) = 0.9453087204829/sqrt(pi);
            case 6
                roots(1) = 0.436077411927617*sqrt(2); weights(1) = 0.7246295952244/sqrt(pi);
                roots(2) = -roots(1);         weights(2) = weights(1);
                roots(3) = 1.335849074013697*sqrt(2); weights(3) = 0.1570673203229/sqrt(pi);
                roots(4) = -roots(3);         weights(4) = weights(3); 
                roots(5) = 2.350604973674492*sqrt(2); weights(5) = 0.00453000995509/sqrt(pi);
                roots(6) = -roots(5);         weights(6) = weights(5);
                
            case 7
                roots(1) = 0.816287882858965*sqrt(2); weights(1) = 0.4256072526101/sqrt(pi);
                roots(2) = -roots(1);         weights(2) = weights(1);
                roots(3) = 1.673551628767471*sqrt(2); weights(3) = 0.05451558281913/sqrt(pi);
                roots(4) = -roots(3);         weights(4) = weights(3); 
                roots(5) = 2.651961356835233*sqrt(2); weights(5) = 0.0009717812450995/sqrt(pi);
                roots(6) = -roots(3);         weights(6) = weights(5);
                roots(7) = 0;                 weights(7) = 0.8102646175568/sqrt(pi);
                
        end
end