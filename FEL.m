%%% FE-ANALYSIS
function [U]=FEL(nelx,nely,x,penal,KE,E,dof)
    K = sparse(2*(nelx+1)*(nely+1),2*(nelx+1)*(nely+1));
    F = sparse(2*(nely+1)*(nelx+1),1); 
    U = zeros(2*(nely+1)*(nelx+1),1);
    for elx = 1:nelx
      for ely = 1:nely
        n1 = (nely+1)*(elx-1)+ely; 
        n2 = (nely+1)* elx   +ely;
        edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
        K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*E(ely,elx)*KE;
      end
    end
    % DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
    F(dof,1) = 1;
    fixeddofs   = union(1:2*(nely+1):(nelx/2)*(nely+1)*2+1,2:2*(nely+1):(nelx/2)*(nely+1)*2+2);
    alldofs     = 1:2*(nely+1)*(nelx+1);
    freedofs    = setdiff(alldofs,fixeddofs);
    % SOLVING
    U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);
    U(fixeddofs,:) = 0;
end