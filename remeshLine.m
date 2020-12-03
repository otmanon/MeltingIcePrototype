%Remeshes domain, making sure to keep boundary of interior F and exterior F intact.
% Input:
% V: Vertices  of geometry
% intF: Faces belonging to solid domain
% extF : Faces belonging to liquid domain 
%Output:
%V2, F2 - new vertices and geometry
% intF - new interior Faces
% extF - new exterior Faces

function [V2, F2, intF2, extF2, corners] = remeshLine(V, intF, extF, maxArea, bboxSize, bboxW, res)
    S1 = boundary_faces(intF);
    S2 = boundary_faces(extF);
    S2Flipped = fliplr(S2);
    interfaceEdges = intersect(S1, S2Flipped, 'rows');  %edges of interface, without the ones in the bounding box
    [RV, g2l, l2g, RS] = remove_unreferenced(V, interfaceEdges);
    [QV, QS] = maintainBoundaryQuality(RV, RS, bboxSize, res);
    [BBV, BBE, junk, corners] = embedLineInBoundingBox(QV, QS, bboxSize, bboxW);
    H = [];
    warning('off')
    [V2,F2,nbs] = triangle(BBV, BBE, H,'MaxArea', maxArea, 'Quality');
    [intF2, extF2] = segmentDomain(V, S1, V2, F2);
end


