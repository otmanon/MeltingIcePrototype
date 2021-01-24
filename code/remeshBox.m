%Remeshes domain, making sure to keep boundary of interior F and exterior F intact.
% Input:
% V: Vertices  of geometry
% intF: Faces belonging to solid domain
% extF : Faces belonging to liquid domain 
%Output:
%V2, F2 - new vertices and geometry
% intF - new interior Faces
% extF - new exterior Faces

function [V2, F2, intF2, extF2, data] = remeshBox(V, intF, extF, maxArea, bboxL, bboxW, avgLength, data, step)
    S1 = boundary_faces(intF);
    S2 = boundary_faces(extF);
    [RV, g2l, l2g, RS] = remove_unreferenced(V, S1);
    [QV, QS] = maintainBoundaryQuality(RV, RS, bboxL, avgLength);
    [extV, extS] = makeBox(bboxL, bboxW, 1);
    Vembedded = [QV; extV];
    extS = extS + size(QS,1 );
    S = [QS; extS];
    H = [];
    warning('off')
    tic;
    [V2,F2,nbs] = triangle(Vembedded, S, H,'MaxArea', maxArea, 'Quality');
    triangle_time = toc
    tic;
    [intF2, extF2] = segmentDomain(V, S1, V2, F2);
    segment_time = toc
    data(step, 3:4) = [triangle_time, segment_time];
end