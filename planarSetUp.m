%sets up planar configuration for classic Stefan Problem/ Mullins Sekerka
%INstability
%Inputs:
%Res: resolution of lines
%MaxArea: Max triangle area
%Outputs:
%V - Vertices of triangulation
%F - Faces of triangulation
%S - interface Edges.
function [V, F, intF, extF] = planarSetUp(res, maxArea, bboxSize)
    startPoint = [0.1 bboxSize];
    stopPoint = [0.1 0];

    [V, E] = buildStraightLine(res, startPoint, stopPoint);
    [Vembedded, Ei, S] = embedInBoundingBox(V, E, bboxSize);
    H = [];
    [V,F,nbs] = triangle(Vembedded, Ei, H,'MaxArea', maxArea, 'Quality');
    [extF, intF ] = segmentDomain(Vembedded, S, V, F);
end

