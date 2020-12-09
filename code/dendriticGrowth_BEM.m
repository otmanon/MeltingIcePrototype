clear; 
addpath('function');
res = 20;       %resolution of line we wish to make. 
maxArea = 1;  %maximum triangle area for remesher
bboxL = 30; %bounding box length
bboxW = 30; %bounding box width

c = 0.05; %surface tension
dt = 0.01; %timestep size
bbT = -10; %boundary temperatue... switch to positive if you wanna see melting
meshScale = 10; % how big should our mesh be.  Usually from 1 to 10 is good depending on size of bbox and input mesh and maxArea
meshName = "../data/pentagonHighRes.obj";

nv = 200;
[gridV,gridF]=create_regular_grid(nv,nv);
gridV = 30*gridV;

[V,S, intS, extS, T, dTdn, avgL] = init(maxArea, bboxL, bboxW, c, bbT, meshScale, meshName, gridV);


minX = 0;
maxX = bboxW;
minY = 0; 
maxY = bboxL; 

map = colormapSetup();
colormap(map);

tsurf(gridF,gridV,'CData',T,fphong,falpha(1,0), 'LineStyle', 'none');
hold on;
colormap(map);
p = plot_edges(V,S,'k','LineWidth',2);
axis([minX maxX minY maxY])

drawnow;



for step = 1:10000
   waitforbuttonpress;
    [V, S, intS, extS, T]=stepSim(V, S, intS, extS, T, dTdn, c, dt, bbT, gridV, avgL);
   
    plot_edges(V,S,'k','LineWidth',2);
    drawnow;

    if mod(step, 10) == 1
        figgif(strcat(strcat("snowflakest", num2str(c)), ".gif"));
    end
    
end

function [V, S, intS, extS, T, dTdn, avgL]=init(maxArea, bboxL, bboxW, c, bbT, scale, meshName, gridV)
    [V, S, intS, extS] = bemSetUp(bboxL, bboxW, scale, meshName);
    avgL = norm(edge_lengths(V, intS));
    [T, dTdn] = temperatureDiffusionBEM(V, S, intS, extS, c, bbT, gridV);
end

function [V2, S, intS, extS, T2]=stepSim(V, S, intS, extS, T, Vn, c, dt, bbT, gridV, avgL)
    b = unique(intS);
    N = per_vertex_normals2D(V, intS);
    V2 = V(b) + Vn.*N*dt;
    V2 = [V2; V(end - 4: end, :)];

    %[V2, F2, intF2, extF2] = remeshBox(V, intF, extF, maxArea, bboxL,bboxW, avgL);
    T2 = temperatureDiffusionBEM(V2, S, intS, extS, c, bbT, gridV);
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


