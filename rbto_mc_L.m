%%% use Monte Carlo
function rbto_mc_L(nelx, nely, penal, x, dismax, dof, a,b)
 
    rng(0);
    nKL = 2;
 
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
 
    numSamples = 100000/2;
    u(1:numSamples) = 0;
    [eigV, eigF] = KL(nelx, nely, nKL);
  
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
       
    data = randn(numSamples,nKL);
%     data = rand(numSamples,nKL);
%     for i = 1:nKL
%         index = randperm(numSamples);
%         prob = (index'-data(:,i))/numSamples;
%         data(:,i) = sqrt(2)*erfinv(2*prob-1);
%     end
    
    tic
    parfor i = 1:numSamples
        Z = sqrt(eigV) .* data(i, :)';
        E(1:nely, 1:nelx) = 0;
        for j = 1:nKL
            E = E + Z(j) * squeeze(eigF(j, :, :));
        end
        E = a + (b - a) * normcdf(E);
        [U] = FEL(nelx, nely, x, penal, KE, E, dof);
        u(i) = U(dof);
    end    
    toc
        
    disp('Prob: '); 
    prob = sum(abs(u) - dismax >= 0)/numSamples
%     disp('1 - Prob: ');
%     1 - sum(dismax - abs(u) <= 0)/numSamples
 
    disp('mean: '); mcsu = mean(u)
    disp('std: '); mcsstd = std(u)
 
    s = length(colPoints);
    v(1:s) = 0;
 
    for i = 1:s
        Z = sqrt(eigV) .* colPoints(i, :)';
        E(1:nely, 1:nelx) = 0;
        for j = 1:nKL
            E = E + Z(j) * squeeze(eigF(j, :, :));
        end
        E = a + (b - a) * normcdf(E);
        [U] = FEL(nelx, nely, x, penal, KE, E, dof);
        v(i) = U(dof);
    end
 
    numVar = 10;
    N(1:s, 1:numVar) = 0;
    for i = 1:s
        N(i, 1) = 1;
        N(i, 2) = colPoints(i, 1);
        N(i, 3) = colPoints(i, 2);
        N(i, 4) = colPoints(i, 1) ^ 2 - 1;
        N(i, 5) = colPoints(i, 2) ^ 2 - 1;
        N(i, 6) = colPoints(i, 1) * colPoints(i, 2);
        N(i, 7) = colPoints(i, 1) ^ 3 - 3 * colPoints(i, 1);
        N(i, 8) = colPoints(i, 2) ^ 3 - 3 * colPoints(i, 2);
        N(i, 9) = colPoints(i, 1) * colPoints(i, 2) ^ 2 - colPoints(i, 1);
        N(i, 10) = colPoints(i, 2) * colPoints(i, 1) ^ 2 - colPoints(i, 2);
    end
 
    a = N'*N \ N' * v';
 
    disp1(1:numSamples) = 0;
    tic
    for i = 1:numSamples
        x = data(i, :);
        disp1(i) = a(1) + a(2) * x(1) + a(3) * x(2) + ...
          a(4) * (x(1) ^ 2 - 1) + a(5) * (x(2) ^ 2 - 1) + ...
          a(6) * x(1) * x(2) + ...
          a(7) * (x(1) ^ 3 - 3 * x(1)) + a(8) * (x(2) ^ 3 - 3 * x(2)) + ...
          a(9) * (x(1) * (x(2) ^ 2) - x(1)) + a(10) * ((x(1) ^ 2) * x(2) - x(2));
    end
    toc
    mdis = mean(disp1)
    stddis = std(disp1)
    
    fileID = fopen('data.txt','a+');
    fprintf(fileID,'\nprob: %4.7f',prob);
    fprintf(fileID,'\nMCS mean: %4.7f',mcsu);
    fprintf(fileID,'\nMCS std: %4.7f',mcsstd);
    fprintf(fileID,'\nSRSM mean: %4.7f',mdis);
    fprintf(fileID,'\nSRSM std: %4.7f',stddis);
    fprintf(fileID,'\n\n');
    fclose(fileID);
    %colormap(gray); imagesc(-x); axis equal; axis tight; axis off; pause(1e-6);

    [f1,x1] = ecdf(u);
    [f2,x2] = ecdf(disp1);
    
    figure
    hold on
    r = length(x1)-10:length(x1);
    plot(x1(r),f1(r),'-ob')
    plot(x2(r),f2(r),'--r*')
    legend('MCS','SRSM','Location','best')
    xlabel('Displacement')
    ylabel('Cumulative probability') 
    box on
    hold off
    
    figure
    hold on
    plot(x1,f1,'-b')
    plot(x2,f2,'--r')
    legend('MCS','SRSM','Location','best')
    xlabel('Displacement')
    ylabel('Cumulative probability')
    box on
    hold off
 
end