% Takes Vertices and edges and finds a surface bounding box to embed them in
function [outV, outE, surfaceE, cornerInds] = embedLineInBoundingBox(V, E, h, w)
    [top, topI] = max(V(:, 2));
    [bottom ,bottomI] = min(V(:, 2));
    n = length(V);
    boundingBox = [0 h; w h; w 0; 0 0];
    newE = [1+n topI; topI 2+n; 2+n 3+n; 3+n bottomI; bottomI 4+n; 4+n 1+n];
    outV = [V; boundingBox];
    outE = [E; newE];
    newBBE = [1+n topI; bottomI 4+n; 4+n 1+n]; %Edges that correspond to boundary of ice (ie square including left edges and ice boundary
    surfaceE = [E; newBBE];
    cornerInds = unique(newBBE);
    
end