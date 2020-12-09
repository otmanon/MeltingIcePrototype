%diffuses temperature, assuming boundary liquid vertices are fixed at
%10degrees, interfaceTemperature is fixed to bulk Melting temperature,
%Using the poisson steady state approximation.

%Inputs:
%V - Defines domain vertices
%S - List of edges that build interior boundary (interface) and
%exterior boundary)
%intS, extS - List of edges that build interior boundary (interface) and
%exterior boundary), separated
%domain, respectively

%Output - Steady state temperature field
% 
function [T, dTdn] = temperatureDiffusionBox(V, S, intS, extS, c, bbT, gridV)
    [Tb, bi] = getBoundaryTempAndIndicesBEM(intS, extS, V, c, bbT);
    [T, dTdn] = poissonSolve(V,S, intS, extS, Tb, bi, gridV);
   % T = projectTempOnBoundary(V, T);
end

function [Tb, bi] = getBoundaryTempAndIndicesBEM(intS, extS, V, surfaceTensionConstant, bbT)
    stateVec = zeros(length(V), 1); %if 0, etnry is interior, if 1, entry is boundary of interior, if 2, entry is boundary of exterior

    exteriorBoundaryIndeces = unique(extS);%(V(:, 1) == rightMostX);
    %exteriorBoundaryIndeces = unique(boundary_faces(extF)); % only set rightmost boundary to temp
    stateVec(exteriorBoundaryIndeces) = 2;

    [interiorTemp, interiorBoundaryIndeces] = gibbsThomsonBoundaryConditions(V, intS, surfaceTensionConstant);
    stateVec(interiorBoundaryIndeces) = 1;
    interiorIndices = find(stateVec == 1);
    exteriorIndices = find(stateVec == 2);
    exteriorTemp = ones(length(exteriorIndices), 1)*bbT;
    Tb = [ interiorTemp; exteriorTemp];
    bi = [  interiorIndices; exteriorIndices];
   
end

function [T, dTdn] = poissonSolve(V, S, intS, extS, Tb, bi, gridV)
   [T, dTdn]=  solveBEM(V, S, intS, extS, Tb, Tb, gridV);
   
end

function [Tb, bi] = gibbsThomsonBoundaryConditions(V, intS, c)
   % cornerInds = unique(corners);
     
    [RV, g2l, l2g, RS] = remove_unreferenced(V, intS);
    
    [k, alpha, ev] = curvature(RV, RS);
   
    Tb = -k * c;
    bi = l2g;
    
end



