%% KL EXPANSION
function [eigV,eigF] = KL(nelx,nely,nKL)

    lenX = nelx+1+nelx*3;
    x = linspace(0,nelx,nelx+1+nelx*3);
    corrLenX = 0.6;
    Cx(1:lenX, 1:lenX) = 0;

    for k = 1:lenX
        for l = 1:lenX
            Cx(k,l) = exp(-abs(x(k) - x(l))/corrLenX);
        end
    end
    
    [lambdax,phix] = IntEqSolver1d(x(1),x(end),lenX,Cx,nKL,'collocation',1,0);
    
    lenY = nely+1+nely*3;
    y = linspace(0,nely,lenY);
    corrLenY = 0.6;
    Cy(1:lenY, 1:lenY) = 0;

    for k = 1:lenY
        for l = 1:lenY
            Cy(k,l) = exp(-abs(y(k) - y(l))/corrLenY);
        end
    end
    
    [lambday,phiy] = IntEqSolver1d(y(1),y(end),lenY,Cy,nKL,'collocation',1,0);
    
    eigV = lambdax.*lambday;
    eigF(1:nKL,1:nely,1:nelx) = 0;
    
    for k = 1:nKL
        xx = phix(3:4:end-2,k);
        yy = phiy(3:4:end-2,k);
        [XX,YY] = meshgrid(xx,yy);
        eigF(k,:,:) = XX.*YY;
    end
    
end