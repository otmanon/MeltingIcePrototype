clear; 

res = 50;       %resolution of line we wish to make. 
maxArea = 0.1;  %maximum triangle area for remesher
bboxH = 10; %bounding box length
bboxW = 20; %bounding box width

c = 0.0; %surface tension
dt = 0.01; %timestep size
bbT = -10; %boundary temperature
numBumps = 1 ;  %Controls number of bumps on the line... keep this a multiple of 2. If you want to make it an odd number, it'll look weird unless you go to disturbedPlanerSetUp.m, and change the cos to a sin in the "perturbStraightLine" function
bumpSize = 0.5; %Controls size of bumps. Keep between 0 to 0.5... otherwise it's too big
[V, F, intF, extF, T, avgL] = init(maxArea, bboxH, bboxW,res, c, bbT, numBumps, bumpSize);


minX = 0;
maxX = bboxW;
minY = 0; 
maxY = bboxH; 

map = colormapSetup();
colormap(map);

p1 = tsurf(F, V, 'CData', T, 'FaceColor','interp',  'LineStyle', 'none');
hold on;
colormap(map);
p2 = tsurf(intF, V, 'CData', T, 'FaceColor','interp');% 'LineStyle', 'none');
axis([minX maxX minY maxY])

drawnow;


waitforbuttonpress;
for step = 1:10000
   
    [V, F, intF, extF, T]=stepSim(V, F, intF, extF, T, c,dt, bbT, maxArea, bboxH, bboxW, avgL);
   
     p1.CData = T;
     p1.Vertices = V;
     p1.Faces = F;
     p2.CData = T;
     p2.Vertices = V;
     p2.Faces = intF;
    drawnow;

    if mod(step, 10) == 1
   %     figgif(strcat(strcat("../gifs/snowflakest", num2str(c)), "numBumps",num2str(numBumps), ".gif"));
    end
    
end

function [V, F, intF, extF, T, avgL]=init(maxArea, bboxL, bboxW, res, c, bbT, numBumps, bumpSize)
    [V, F, intF, extF, avgL] = disturbedPlanarSetUp(maxArea, bboxL, bboxW, res, numBumps, bumpSize);

    T = temperatureDiffusionLine(V, F, intF, extF, c, bbT);
end

function [V2, F2, intF2, extF2, T2]=stepSim(V, F, intF, extF, T, c, dt, bbT, maxArea, bboxL, bboxW, avgL)
    TGrad = getTemperatureGradient(V, F, T);
    
    [dTdn, N, M, S] = getFluxAlongBoundary(V,intF, F, TGrad);
    VMotion = fitVertexMotion(V, S, N, dTdn);
    V = V + VMotion*dt;
   

    [V2, F2, intF2, extF2] = remeshLine(V, intF, extF, maxArea, bboxL,bboxW, avgL);
%     V2 = V;
%     F2 = F;
%     intF2 = intF;
%     extF2 = extF;

    T2 = temperatureDiffusionLine(V2, F2, intF2, extF2, c, bbT);
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


