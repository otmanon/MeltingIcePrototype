%Takes vertices with edges connecting them, ensures the edges are within a
%range of distances with each other, otherwise subdivides
%   Detailed explanation goes here
function [QV, QS] = maintainBoundaryQuality(V, F,  avgLength)
    maxLength = avgLength + avgLength/2;
    minLength = avgLength - avgLength/1.1;
    [V1, F1] = subdivideLargeEdges(V, F, maxLength);
   % [V2, F2] = mergeSmallEdges(V1, F, minLength);
%=      [V3, S3] = subdivideLargeEdges(V2, S2, maxLength);
%      [V4, S4] = mergeSmallEdges(V3, S3, minLength);
    QV = V1;
    QS = F1;
end

function [Vout, Fout] = subdivideLargeEdges(V, F, maxLength)
    S = edges(F);
    Vout = V;
    Fout = F;
    D =  V(S(:, 2), :) - V(S(:, 1), :);
    D = vecnorm(D, 2, 2);
    largeInd = find(D > maxLength); %Contains idnex into S of edges that must be subdivided
    
    lE = S(largeInd, :);
    if (size(lE, 1) > 0)
        [Vout, Fout, bF] = split_edges(V, F, lE);
    end;
%     numOldV = size(V, 1);
%     numNewV = size(largeInd, 1);
%     if (size(largeInd, 1) > 0)
%         newEdgeMidpoints = (V(S(largeInd(:), 2), :) + V(S(largeInd(:), 1), :))/2;
%     
%         newVIndices = numOldV+1:numOldV+numNewV;
% 
%         oldEdgeCopy = S(largeInd, :);
% 
%         S(largeInd, 1) = newVIndices;
%         oldEdgeCopy(:, 2) = newVIndices;
% 
%         Sout = [S; oldEdgeCopy];
%         Vout = [V; newEdgeMidpoints];
%     else
%         Sout = S;
%         Vout = V;
%     end
% 
%     
end


function [Vout, Fout] = mergeSmallEdges(V, F, minLength)
    [Vout, Fout] = collapse_close_points(V, F, minLength);
  
   
end

