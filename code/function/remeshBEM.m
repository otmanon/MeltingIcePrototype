%Remeshes domain, making sure to keep boundary of interior F and exterior F intact.
% Input:
% V: Vertices  of geometry
% intF: Faces belonging to solid domain
% extF : Faces belonging to liquid domain 
%Output:
%V2, F2 - new vertices and geometry
% intF - new interior Faces
% extF - new exterior Faces

function [V2, S2, intS2, extS2] = remeshBox(V,S, intS, extS, bboxL, bboxW, avgLength)
    

    [QV, QS] = maintainBoundaryBEM(V(1:end-4, :), intS, bboxL, avgLength);
    extS2 = extS + size(QV, 1) - size(V, 1) + 4;
    QV = [QV; V(end-3:end, :)];
    V2 = QV;

    intS2 = QS;
        
    S2 = [intS2; extS2];
    

end