function [lambda,phi] = IntEqSolver1d(a,b,n,Cd,neig,method,p,ToPlot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% returns solution of the Fredholm integral equation of the second kind
% [a,b] - interval of integration
% n     - grid size = total number of nodes including endpoints
% neig   - number of eigenpairs required
%
% method: 'collocation' or 'galerkin'
% p - degree of the polynomials in Lagrange basis: either 1 or 2
%
% varargin - optional list of arguments 
% contains parameters for the kernels that need it: eta, sigma
%
% examples of run:
% [lambda,phi] = IntEqSolver1d(0,1,50,10,'exponential','collocation',2,1,1/10,1)
% Written by Veronika Vasylkivska <vasylkiv@math.oregonstate.edu>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off all

h = (b-a)/(n-1);

if p == 1          %% piecewise linear functions are used for basis
    nodes = a:h:b; %% grid nodes
    nnodes = n;
       
elseif p == 2        %% piecewise quadratic functions are used for basis
    nodes = a:h/2:b; %% grid nodes
    nnodes = 2*n-1;
else
    error('Not valid degree of polynomials: should be 1 or 2');
end

% find matrices for the generalized eigenvalue problem
% using one of the methods: collocation or Galerkin

D = zeros(nnodes,nnodes);
switch lower(method)
    case 'collocation'
        %% set up number of integration points nw, nodes xw, and weights w
        nw = 11;
        [xw,w] = GaussQuad(nw,'legendre');
                
        %% for the first element the weights and roots are different
        xw1 = 0.5*(xw+1); w1 = 0.5*w;
        C = StatCov1d(nodes,h*xw1+nodes(1),nodes,Cd,0);                
        psi = shapefun(xw1,p,1); %% calculations on ref.element
        D(:,1) = C*(w1.*psi);

        %% for the last element the weights and roots are different as well
        xw2 = 0.5*(xw-1); w2 = w1;
        C = StatCov1d(nodes,h*xw2+nodes(nnodes),nodes,Cd,0);                
        psi = shapefun(xw2,p,1); %% calculations on ref.element
        D(:,nnodes) = C*(w2.*psi);

        psi = shapefun(xw,p,1); %% calculations on ref.element
        if p == 2; psi2 = shapefun(xw,p,0); end

        for j = 2:nnodes-1
            C = StatCov1d(nodes,h*xw+nodes(j),nodes,Cd,0);
            if p == 1
                D(:,j) = C*(w.*psi);
            else
                if mod(j,2)==0
                    D(:,j) = C*(w.*psi2);
                else
                    D(:,j) = C*(w.*psi);
                end
            end
        end
                
        D = 2*h*D;
        L = eye(nnodes,nnodes);
                  
    case 'galerkin'
        nw = 11;
        [xw,w] = GaussQuad(nw,'legendre');
        xw1 = 0.5*(xw+1); w1 = 0.5*w; %% weights and roots for the first element
        xw2 = 0.5*(xw-1); w2 = w1; %% weights and roots for the last element
        
        psi1 = shapefun(xw1,p,1); psi2 = shapefun(xw2,p,1); %% shape function on the 1st and last element
        
        %% shape functions on other elements
        if p == 1
            psi(:,1) = shapefun(xw,p,1);
            psi(:,2) = psi(:,1);
        elseif p == 2
            psi(:,1) = shapefun(xw,p,0); psi(:,2) = shapefun(xw,p,1); 
        end       
        
        %% calculations of the matrix D entries
        W = w1*w1'; 
        C = StatCov1d(h*xw1+nodes(1),h*xw1+nodes(1),nodes,Cd,0);
        D(1,1) = sum(sum(C.*W.*(psi1*psi1')));

        W = w1*w2';
        C = StatCov1d(h*xw1+nodes(1),h*xw2+nodes(nnodes),nodes,Cd,0);         
        D(1,nnodes) = sum(sum(C.*W.*(psi1*psi2')));
        D(nnodes,1) = D(1,nnodes);

        W = w1*w';
        for j = 2:nnodes-1
            C = StatCov1d(h*xw1+nodes(1),h*xw+nodes(j),nodes,Cd,0);
            D(1,j) = sum(sum(C.*W.*(psi1*psi(:,mod(j,2)+1)')));
            D(j,1) = D(1,j);
        end

        W = w2*w2';
        C = StatCov1d(h*xw2+nodes(nnodes),h*xw2+nodes(nnodes),nodes,Cd,0);         
        D(nnodes,nnodes) = sum(sum(C.*W.*(psi2*psi2')));
            
        W = w2*w';
        for j = 2:nnodes-1
            C = StatCov1d(h*xw2+nodes(nnodes),h*xw+nodes(j),nodes,Cd,0);
            D(nnodes,j) = sum(sum(C.*W.*(psi2*psi(:,mod(j,2)+1)')));
            D(j,nnodes) = D(nnodes,j);
        end    
            
        W = w*w';
        for j = 2:nnodes-1
            for k = j:nnodes-1
                C = StatCov1d(h*xw+nodes(j),h*xw+nodes(k),nodes,Cd,0);
                D(j,k) = sum(sum(C.*W.*(psi(:,mod(j,2)+1)*psi(:,mod(k,2)+1)')));
                D(k,j) = D(j,k);
            end
        end
        D = 4*h*h*D;
        L = GalerkinL(n,h,p);
        
end

% get eigenpairs
[phi,lambda] = eigs(D,L,neig);

% improve approximation of the eigenfunctions
res = -(D*phi-L*phi*lambda);  %% find residual first
lambda = diag(lambda,0);

% for j = 1:neig
%    M = D - lambda(j)*L;
%    phi(:,j) = phi(:,j) + M\res(:,j);
% end

% get rid of 0 complex part (needed for gaussian kernels)
lambda = abs(real(lambda));

% sort found eigenvalues in descending order
[lambda,ix] = sort(lambda,'descend');

% rearrange the eigenfunctions in the corresponding order
Nphi = phi;            
for j = 1:neig
    phi(:,ix(j)) = Nphi(:,j);
end 
phi = signchange(phi);

% normalize found eigenfunctions
if p == 1
    for j = 1:neig
        v = phi(:,j);
        vnorm = sqrt(h/3*(v(1)^2+v(n)^2+2*sum(v(2:n-1).*v(2:n-1))+sum(v(1:n-1).*v(2:n))));
        phi(:,j) = phi(:,j)/vnorm;
    end
else
    for j = 1:neig
        v = phi(:,j);
        vnorm = sqrt(h/15*(2*v(1)^2+2*v(nnodes)^2+8*sum(v(2:2:nnodes-1).^2)+4*sum(v(3:2:nnodes-1).^2)+...
            2*sum(v(1:nnodes-1).*v(2:nnodes))-sum(v(1:2:nnodes-2).*v(3:2:nnodes))));
        phi(:,j) = phi(:,j)/vnorm;
    end
end

% make sure to get rid of complex zero part for some kernels
for j = 1:neig
    phi(:,j) = real(phi(:,j));
end

% plot eigenvalues if needed
if ToPlot
    plot(lambda,'-sk','LineWidth',1.5);
end

end
       
%%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = shapefun(x,p,mode) %%%%% shape function for reference element

    lenx = length(x);
    f = zeros(lenx,1);
    
    if p == 1 %% piecewise linear
        for j = 1:lenx
            if (x(j)<0)&&(x(j)>=-1)
                f(j) = 1+x(j);
            elseif (x(j)>=0)&&(x(j)<=1)
                f(j) = 1-x(j);
            else
                f(j) = 0;
            end
        end
        
    elseif p == 2 %% piecewise quadratic
        switch mode
        case 1 %% mesh vertices
            for j = 1:lenx
                if (x(j)<0)&&(x(j)>=-1)
                    f(j) = (1+x(j)).*(1+2*x(j));
                elseif (x(j)>=0)&&(x(j)<=1)
                    f(j) = (1-x(j)).*(1-2*x(j));
                else
                    f(j) = 0;
                end
            end
        case 0 %% midpoints
            for j = 1:lenx
                if abs(x(j)) < 0.5
                    f(j) = 1-4*x(j).*x(j);
                else
                    f(j) = 0;
                end
            end    
        end
    end
end


function ResF = signchange(F)
    ResF = F;
    n = size(F,2);
    for j = 1:n
        if F(1,j) < 0
            ResF(:,j) = -F(:,j);
        end
    end
end


function L = GalerkinL(n,h,p)

    switch p
        case 1
            L = zeros(n,n);
            
            
            L(1,1) = 1/3;
            L(n,n) = 1/3;
            for k = 2:n-1
                L(k,k) = 2/3;
            end
            for k = 1:n-1
                L(k,k+1) = 1/6;
                L(k+1,k) = 1/6;
            end
            L = h*L;
%             A = h/6*[2 1; 1 2];
%             for k = 1:n-1
%                 for i = 1:2
%                     for j = 1:2
%                         ig = k+i-1;
%                         jg = k+j-1;
%                         L(ig,jg) = L(ig,jg) + A(i,j);
%                     end
%                 end
%             end
        case 2
            L = zeros(2*n-1,2*n-1);
            A = h/30*[4 2 -1; 2 16 2; -1 2 4];
            
            for k = 1:n-1
                for i = 1:3
                    for j = 1:3
                        ig = 2*k+i-2;
                        jg = 2*k+j-2;
                        L(ig,jg) = L(ig,jg) + A(i,j);
                    end
                end
            end
    end    
end

function M = StatCov1d(x,y,xdata,C,ToPlot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% returns the covariance function 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y] = meshgrid(xdata);

M = interp2(X,Y,C,x,y);
M = M';

% plot the covariance if needed
if ToPlot
    figure;
    surf(x,y,M);
    title(strcat(kernel,' model'),'fontsize',20);
end
end