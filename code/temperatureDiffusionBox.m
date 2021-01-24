%diffuses temperature, assuming boundary liquid vertices are fixed at
%10degrees, interfaceTemperature is fixed to bulk Melting temperature,
%Using the poisson steady state approximation.

%Inputs:
%V, F - Defines domain geometry and topolgy
%S - boundary of Solid domain
%intF, extF - List of faces that are inside solid domain and outside solid
%domain, respectively

%Output - Steady state temperature field
% 
function T = temperatureDiffusionBox(V, F, intF, extF, c, bbT)
    [Tb, bi] = getBoundaryTempAndIndicesBox(intF, extF, V, c, bbT);
    T = poissonSolve(V, F, Tb, bi);
   % T = projectTempOnBoundary(V, T);
end

function [Tb, bi] = getBoundaryTempAndIndicesBox(intF, extF, V, surfaceTensionConstant, bbT)
    stateVec = zeros(length(V), 1); %if 0, etnry is interior, if 1, entry is boundary of interior, if 2, entry is boundary of exterior
    rightMostX = max(V(:, 1));
    exteriorBoundaryIndeces = unique(boundary_faces(extF));%(V(:, 1) == rightMostX);
    %exteriorBoundaryIndeces = unique(boundary_faces(extF)); % only set rightmost boundary to temp
    stateVec(exteriorBoundaryIndeces) = 2;

    [interiorTemp, interiorBoundaryIndeces] = gibbsThomsonBoundaryConditions(V, intF, surfaceTensionConstant);
    stateVec(interiorBoundaryIndeces) = 1;
    interiorIndices = find(stateVec == 1);
    exteriorIndices = find(stateVec == 2);
    exteriorTemp = ones(length(exteriorIndices), 1)*bbT;
    Tb = [ interiorTemp; exteriorTemp];
    bi = [  interiorIndices; exteriorIndices];
   
end

function T = poissonSolve(V, F, Tb, bi)
    C = cotmatrix(V, F); %Quadratic Coefficients
    B = []; %linear Coefficients empty
    [T, prefactoredSystem] = min_quad_with_fixed(C, B, bi, Tb);
   
end

function [Tb, bi] = gibbsThomsonBoundaryConditions(V, intF, c)
   % cornerInds = unique(corners);
    mode =4;
   
    S = boundary_faces(intF);   

    [RV, g2l, l2g, RS] = remove_unreferenced(V, S);
         BC = barycenter(RV, RS);
    [k, alpha, ev] = curvature(RV, RS);
    eNormals = perEdgeNormals(RV, RS);
    eAxis = repmat([1, 0], size(eNormals, 1), 1);
    theta = acos(dot(eAxis', eNormals'));
    sigma = c.*( 1 - cos((theta).*mode));
    sigmaV = zeros(size(sigma, 1), 1);
    sigmaV(RS(:, 1)) = sigma;
    sigmaV(RS(:, 2)) = sigmaV(RS(:, 2)) + sigma;
    Tb = -k .* sigmaV';
    n = eNormals;
    n = n .* c; %sigma';
    
    %quiver(BC(:, 1), BC(:, 2), n(:, 1), n(:, 2));
    bi = l2g;
    
end


function T = projectTempOnBoundary(V,  T)
    eps = 0.1;
    topY = max(V(:, 2)); 
    botY = min(V(:, 2)); 
    leftX = min(V(:, 1));
    
    i1 = find((V(:, 2) < topY + eps) & (V(:, 2) > topY - eps));
    i2 = find((V(:, 2) < botY + eps) & (V(:, 2) > botY - eps));
    i3 = find((V(:, 1) < leftX + eps) & (V(:, 1) > leftX - eps));
    ind = unique([i1; i2; i3]);
    T(ind) = 0;
end
