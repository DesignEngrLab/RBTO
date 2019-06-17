%%% dto(60,20,3,1.5,dismax)
function dto(nelx,nely,penal,rmin,dismax)
    %%% Initial values 
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
 
    dof = 2;
    
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
      
    fval(1,1) = 1;
    change = 1.;
    loop = 0;
    E(1:nely,1:nelx) = 1.35;
    while (change > 0.001)
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
    
%     fig = figure;
%     colormap(gray);
%     imagesc(1-xphy);
%     axis equal; 
%     axis tight; 
%     axis off;
%     
%     fig.PaperPositionMode = 'auto';
%     fig_pos = fig.PaperPosition;
%     fig.PaperSize = [fig_pos(3) fig_pos(4)];
%     print(fig,'MySavedFile','-dbmp','-r0')
    disp('END OPTIMIZATION');
    %rbto_mc(nelx, nely, penal, xphy, dismax, dof, 1, 1.5);
end