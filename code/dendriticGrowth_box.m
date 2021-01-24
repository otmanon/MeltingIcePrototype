clear; 
addpath('function');
res = 20;       %resolution of line we wish to make. 
maxArea = 1;  %maximum triangle area for remesher
bboxL = 40; %bounding box length
bboxW = 40; %bounding box width

c = 0.00; %surface tension
dt = 0.01; %timestep size
bbT = -10; %boundary temperatue... switch to positive if you wanna see melting
meshScale = 2; % how big should our mesh be.  Usually from 1 to 10 is good depending on size of bbox and input mesh and maxArea
meshName = "../data/circle.obj"
[V, F, intF, extF, T, avgL] = init(maxArea, bboxL, bboxW, c, bbT, meshScale, meshName);


minX = 0;
maxX = bboxW;
minY = 0; 
maxY = bboxL; 

map = colormapSetup();
colormap(map);

p1 = tsurf(F, V, 'CData', T, 'FaceColor','interp',  'LineStyle', 'none');
hold on;
colormap(map);
p2 = tsurf(intF, V, 'CData', T, 'FaceColor','interp');% 'LineStyle', 'none');
axis([minX maxX minY maxY])

drawnow;

start_time = now;
timing_data = zeros(10000, 4);
waitforbuttonpress;
for step = 1:10000
    tic
    [V, F, intF, extF, T, timing_data]=stepSim(V, F, intF, extF, T, c,dt, bbT, maxArea, bboxL, bboxW, avgL, timing_data, step);
    time_step = toc;
    timing_data(step, 5) = time_step;
    timing_data(step, 6) = length(boundary_faces(intF));
     p1.CData = T;
     p1.Vertices = V;
     p1.Faces = F;
     p2.CData = T;
     p2.Vertices = V;
     p2.Faces = intF;
    drawnow;

    if mod(step, 10) == 1
       % figgif("../gifs/circleIsotropicMode4st005.gif");
    end
    
end

function [V, F, intF, extF, T, avgL]=init(maxArea, bboxL, bboxW, c, bbT, scale, meshName)
    [V, F, intF, extF] = boxSetUp(maxArea, bboxL, bboxW, scale, meshName);
    avgL = mean(edge_lengths(V, boundary_faces(intF)));
    T = temperatureDiffusionBox(V, F, intF, extF, c, bbT);
end

function [V2, F2, intF2, extF2, T2, data]=stepSim(V, F, intF, extF, T, c, dt, bbT, maxArea, bboxL, bboxW, avgL, timing_data, step)
    TGrad = getTemperatureGradient(V, F, T);
    
    [dTdn, N, M, S] = getFluxAlongBoundary(V,intF, F, TGrad);
    VMotion = fitVertexMotion(V, S, N, dTdn);
    V = V + VMotion*dt;
   
    tic;
    [V2, F2, intF2, extF2, data] = remeshBox(V, intF, extF, maxArea, bboxL,bboxW, avgL, timing_data, step);
    time_remesh = toc;
    tic
    T2 = temperatureDiffusionBox(V2, F2, intF2, extF2, c, bbT);
    time_diffusion = toc;
    data = timing_data;
    data(step, 1:2) = [time_remesh time_diffusion];
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


