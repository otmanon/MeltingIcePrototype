clear; 
addpath('function');

maxArea = 1;  %maximum triangle area for remesher
bboxL = 200; %bounding box length
bboxW = 200; %bounding box width

c = 0; %surface tension
dt = 0.01; %timestep size
bbT = -10; %boundary temperatue... switch to positive if you wanna see melting
meshScale = 20; % how big should our mesh be.  Usually from 1 to 10 is good depending on size of bbox and input mesh and maxArea
meshName = "../data/hexagonHighRes.obj";

nv = 1;
[gridV,gridF]=create_regular_grid(nv,nv);
gridV = bboxL*gridV;

[V,S, intS, extS, T, dTdn, avgL] = init(maxArea, bboxL, bboxW, c, bbT, meshScale, meshName, gridV);
minX = 0;
maxX = bboxW;
minY = 0; 
maxY = bboxL; 


map = colormapSetup();
colormap(map);
set(gcf,  'Position', [10, 10, 1200, 1000]);
hold on;
%tsurf(gridF,gridV,'CData',T,fphong,falpha(1,0), 'LineStyle', 'none');

colormap(map);

b = unique(intS);
p = plot_edges(V,S,'k','LineWidth',2);
hold on;
N = perEdgeNormals(V, intS).* dTdn;
BC = barycenter(V, S);
quiver(BC(1:end-4, 1), BC(1:end-4, 2), N(:, 1), N(:, 2)) ;
axis([minX maxX minY maxY])

drawnow;


%waitforbuttonpress;
for step = 1:10000
    cla reset
   % tsurf(gridF,gridV,'CData',T,fphong,falpha(1,0), 'LineStyle', 'none');
    hold on;
    p = plot_edges(V,S,'k','LineWidth',2);
    [V, S, intS, extS, T, dTdn]=stepSim(V, S, intS, extS, T, dTdn, c, dt, bbT, gridV, bboxL, bboxW, avgL);
    
    

    
    
    drawnow;

    if mod(step, 10) == 1
        figgif(strcat(strcat("BEMsnowflake2", num2str(c)), ".gif"));
    end
    
end

function [V, S, intS, extS, T, dTdn, avgL]=init(maxArea, bboxL, bboxW, c, bbT, scale, meshName, gridV)
    [V, S, intS, extS] = bemSetUp(bboxL, bboxW, scale, meshName);
    avgL = (edge_lengths(V, intS));
    avgL = mean(avgL);
    [T, dTdn] = temperatureDiffusionBEM(V, S, intS, extS, c, bbT, gridV);
end

function [V2, S2, intS2, extS2, T2, dTdn2]=stepSim(V, S, intS, extS, T, dTdn, c, dt, bbT, gridV, bboxL, bboxW, avgL)
    b = unique(intS);
    N = perEdgeNormals(V, intS);
    
    VMotion = fitVertexMotion(V, intS, N, dTdn);
    V2 = V(b, :) + VMotion(b, :)*dt;
    V2 = [V2; V(end - 3: end, :)];

    [V2, S2, intS2, extS2] = remeshBEM(V2,S , intS, extS,bboxL,bboxW, avgL);
    [T2, dTdn2] = temperatureDiffusionBEM(V2, S2, intS2, extS2, c, bbT, gridV);
end
function map = colormapSetup()
    teal = [0.3 0.8 0.9];
    orange = [0.8 0.4 0.2];
    deepBlue = [0.3 0.4 1.0];
    colorMapRes = 20;
    rVals = linspace(deepBlue(1), teal(1), colorMapRes);
    gVals = linspace(deepBlue(2), teal(2),  colorMapRes);
    bVals = linspace(deepBlue(3), teal(3), colorMapRes);

    map = [rVals; gVals; bVals]';
end


