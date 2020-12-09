function [V, S, intS, extS] = bemSetUp( bboxL, bboxW, scale, meshName)
    [intV, intF] = readOBJ(meshName);
    V(:, 1) = (intV(:, 1) - min(intV(:, 1)))/(max(intV(:, 1)) - min(intV(:, 1)));
    V(:, 2) = (intV(:, 2) - min(intV(:, 2)))/(max(intV(:, 2)) - min(intV(:, 2)));
    V = V*scale;
    intS = boundary_faces(intF);
    [RV, g2l, l2g, RS] = remove_unreferenced(V, intS);
    RV = RV + [bboxL/2-0.5*scale bboxW/2-0.5*scale];
    [extV, extS] = makeBox(bboxL, bboxW, 1);
    extS = extS + size(RV, 1);
    V = [RV; extV];
    S = [RS; extS];
   	
    
end

