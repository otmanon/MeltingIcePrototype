%Takes vertices with edges connecting them, ensures the edges are within a
%range of distances with each other, otherwise subdivides
%   Detailed explanation goes here
function [QV, QS] = maintainBoundaryBEM(V, S, bbsize, avgLength)
    maxLength = avgLength + avgLength/2;
    minLength = avgLength - avgLength/1.1;
    [V1, S1] = subdivideLargeEdges(V, S, maxLength);
   % [V2, S2] = mergeSmallEdges(V1, S1, minLength);
%=      [V3, S3] = subdivideLargeEdges(V2, S2, maxLength);
%      [V4, S4] = mergeSmallEdges(V3, S3, minLength);
    QV = V1;
    QS = S1;
end

function [Vout, Sout] = subdivideLargeEdges(V, S, maxLength)
    D =  V(S(:, 2), :) - V(S(:, 1), :);
    D = vecnorm(D, 2, 2);
    largeInd = find(D > maxLength); %Contains idnex into S of edges that must be subdivided
    
    numOldV = size(V, 1);
    numNewV = size(largeInd, 1);
    if (size(largeInd, 1) > 0)
        newEdgeMidpoints = (V(S(largeInd(:), 2), :) + V(S(largeInd(:), 1), :))/2;
    
        newVIndices = (numOldV+1):(numOldV+numNewV);

        oldEdgeCopy = S(largeInd, :);

        S(largeInd, 1) = newVIndices;
        oldEdgeCopy(:, 2) = newVIndices;

        Sout = [S; oldEdgeCopy];
        Vout = [V; newEdgeMidpoints];
    else
        Sout = S;
        Vout = V;
    end

    
end


function [Vout, Sout] = mergeSmallEdges(V, S, minLength)
    [Vout, Sout] = collapse_close_points(V, S, minLength);
  
   
end

