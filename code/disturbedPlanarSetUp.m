
function [V, F, intF, extF, avgL] = disturbedPlanarSetUp( maxArea, bboxL, bboxW, res, numBumps, bumpSize)
    startPoint = [bboxW/10 bboxL];
    stopPoint = [bboxW/10 0];

    [V, E] = buildStraightLine(res, startPoint, stopPoint);
    V = perturbStraightLine(V,bboxL, bboxW, numBumps, bumpSize);
    avgL =  norm(edge_lengths(V, E));
    [Vembedded, Ei, S, corners] = embedLineInBoundingBox(V, E, bboxL, bboxW);
    H = [];
    [V,F,nbs] = triangle(Vembedded, Ei, H,'MaxArea', maxArea, 'Quality');
    [extF, intF ] = segmentDomain(Vembedded, S, V, F);
end

%Takes vertices corresponding to a straight line
function V = perturbStraightLine(V,  bboxL,bboxW, numBumps, bumpSize)
    minX = 0;
    maxX = bboxW;
    minY = bboxL/2 - bboxL / 4;
    maxY = bboxL/2 + bboxL / 4;
    
    range = (V(:, 2) < maxY) & (V(:, 2) > minY);
    wavelength = (maxY - minY);
    f = numBumps/ (wavelength);
    amount = 0;
    if (mod(numBumps, 2) == 0)
        amount = -range.*sin(f*V(:, 2)*pi)*bumpSize;
    else
         amount = -range.*cos(f*V(:, 2)*pi)*bumpSize;
    end;
    V(:, 1) = V(:, 1) + amount;
end