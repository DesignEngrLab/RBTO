%%% rbto_den(60,20,3,1.5,dismax)
function rbto_den(nelx,nely,penal,rmin,dismax)

    %%% Initial values
    nKL = 2;
    
    init_val = 0.5;
    x(1:nely,1:nelx) = init_val; 
    xphy(1:nely,1:nelx) = init_val; 

    nu = 0.3;
    k = [ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ... 
         -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];
    KE = 1/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
                      k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                      k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                      k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                      k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                      k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                      k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                      k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
    %%% MMA
    m = 1;
    n = nelx*nely;
    xmin(1:n,1) = 0.001;
    xmax(1:n,1) = 1;
    low = xmin;
    upp = xmax;
    a0 = 1;
    a(1:m,1) = 0;   
    ci(1:m,1) = 1000;
    d(1:m,1) = 0;
    xold1 = x;
    xold2 = x;
    
    [eigV,eigF] = KL(nelx,nely,nKL);
    
    dof = 2;
    upE = 1;
    lwE = 1.5;

    roots = [sqrt(3 + sqrt(6))...
            -sqrt(3 + sqrt(6))...
             sqrt(3 - sqrt(6))...
            -sqrt(3 - sqrt(6))]; 
    colPoints = [ 0 0;...
                  roots(1) 0;...
                  0 roots(1);...
                  roots(2) 0;...
                  0 roots(2);...
                  roots(3) 0;...
                  0 roots(3);...
                  roots(4) 0;...
                  0 roots(4);...
                  roots(3) roots(3);...
                  roots(4) roots(4);...
                  roots(3) roots(4);...
                  roots(4) roots(3);...
                  roots(1) roots(3);...
                  roots(1) roots(4);...
                  roots(3) roots(1);...
                  roots(4) roots(1);...
                  roots(2) roots(3);...
                  roots(3) roots(2);...
                  roots(2) roots(4);...
                  roots(4) roots(2);...
                  roots(1) roots(2);...
                  roots(2) roots(1)];

    %%% Prepare filter
    iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
    jH = ones(size(iH));
    sH = zeros(size(iH));
    k = 0;
    for i1 = 1:nelx
      for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
          for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
            e2 = (i2-1)*nely+j2;
            k = k+1;
            iH(k) = e1;
            jH(k) = e2;
            sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
          end
        end
      end
    end
    H = sparse(iH,jH,sH);
    Hs = sum(H,2);
    
    mpp = [0 0];
    mpptol = 1.;    
    
    fval(1,1) = 1;
    main_loop = 0;
    while mpptol > 0.001
        mppold = mpp;
        main_loop = main_loop + 1;
        loop = 0;
        change = 1.;
        x(1:nely,1:nelx) = init_val;
        xphy = x;
        E = update_E(eigV, eigF, mpp, nKL, upE, lwE, nelx, nely);
        
        while change > 0.001
            loop = loop + 1;
            [U] = FE(nelx,nely,xphy,penal,KE,E,dof);

            dcf(1:nely,1:nelx) = 0; 
            for ely = 1:nely
                for elx = 1:nelx
                      n1 = (nely+1)*(elx-1)+ely; 
                      n2 = (nely+1)* elx   +ely;
                      Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
                      dcf(ely,elx) = -Ue'*penal*xphy(ely,elx)^(penal-1)*E(ely,elx)*KE*Ue;
                end
            end
            c = sum(xphy(:));
            dc(1:nely,1:nelx) = 1; 
            dc(:) = H*(dc(:)./Hs);

            f0val = c;
            df0dx = dc(:);
            fval(1,1) = U(dof)/dismax -1;
            dcf(:) = H*(dcf(:)./Hs);
            dfdx(1,1:n) = dcf(:)/dismax;

            %%% The MMA subproblem is solved at the point xval:
            outeriter = loop;
            [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
            mmasub(m,n,outeriter,x(:),xmin,xmax,xold1(:),xold2(:), ...
            f0val,df0dx,fval,dfdx,low,upp,a0,a,ci,d);

            %%% Update            
            xold2 = xold1;
            xold1 = x;
            xnew = reshape(xmma,nely,nelx);
            xphy(:) = (H*xnew(:))./Hs;
            
            change = max(abs(xnew(:)-x(:)));
            x = xnew;

            % Results
            disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
                  ' Vol Frac.: ' sprintf('%6.4f',sum(sum(xphy))/(nelx*nely)) ...
                  ' ch.: ' sprintf('%6.3f',change )])
            % PLOT DENSITIES  
            colormap(gray); 
            imagesc(1-xphy);
            axis equal; 
            axis tight; 
            axis off;
            pause(1e-6);
        end
        
        mpp = find_mpp(colPoints, eigV, eigF, nKL, nelx, nely, xphy, penal,...
                       KE, dof, upE, lwE, dismax);
        mpptol = max(abs(mpp(:)-mppold(:)));
    end
    
    %%%
    main_loop
    disp('END OPTIMIZATION');
    rbto_mc(nelx, nely, penal, xphy, dismax, dof, upE, lwE);

end

function newE = update_E(eigV, eigF, mpp, nKL, c, d, nelx, nely)

    Z = sqrt(eigV).* mpp';
    E(1:nely,1:nelx) = 0;
    for j = 1:nKL
        E = E + Z(j)*squeeze(eigF(j,:,:));
    end
    newE =  c + (d-c)*normcdf(E);

end

function mpp = find_mpp(colPoints, eigV, eigF, nKL, nelx, nely, x, penal,...
                        KE, dof, c, d, dismax)

    s = length(colPoints);
    v(1:s) = 0;
    
    for i = 1:s
        Z = sqrt(eigV).*colPoints(i,:)';
        E(1:nely,1:nelx) = 0;
        for j = 1:nKL
            E = E + Z(j)*squeeze(eigF(j,:,:));
        end
        E =  c + (d-c)*normcdf(E);
        [U] = FE(nelx,nely,x,penal,KE,E,dof);
        v(i) = U(dof);
    end
    
    numVar = 10;
    N(1:s,1:numVar) = 0;
    for i = 1:s
        N(i,1) = 1;
        N(i,2) = colPoints(i, 1);
        N(i,3) = colPoints(i, 2);
        N(i,4) = colPoints(i, 1)^2 - 1;
        N(i,5) = colPoints(i,2)^2 - 1;
        N(i,6) = colPoints(i,1)*colPoints(i,2);
        N(i,7) = colPoints(i,1)^3 - 3*colPoints(i,1);
        N(i,8) = colPoints(i,2)^3 - 3*colPoints(i,2);
        N(i,9) = colPoints(i,1)*colPoints(i,2)^2 - colPoints(i,1);
        N(i,10) = colPoints(i,2)*colPoints(i,1)^2 - colPoints(i,2);
    end
    
    a = N'*N \ N'*v';
    
    h = @(x) a(1) + a(2)*x(:,1) + a(3)*x(:,2) +...
             a(4)*(x(:,1)^2 - 1) + a(5)*(x(:,2)^2 - 1) +...
             a(6)*x(:,1)*x(:,2) +...
             a(7)*(x(:,1)^3 - 3*x(:,1)) + a(8)*(x(:,2)^3 - 3*x(:,2)) +...
             a(9)*(x(:,1)*(x(:,2)^2) - x(:,1)) + a(10)*((x(:,1)^2)*x(:,2) - x(:,2));
    %g = @(x) h(x)/dismax - 1;
    g = @(x) dismax/h(x) - 1;
    %T = @(x) x;
    Tinv = @(u) u;

    res = CODES.reliability.iform(g,2,2.5,'solver','hmv','Tinv',Tinv);
    disp(res)
    mpp = res.MPTP;

end