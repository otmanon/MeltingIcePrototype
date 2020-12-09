function [N] = per_vertex_normals2D(V, S)
%PER_VERTEX_NORMALS2D Summary of this function goes here
%   Detailed explanation goes here
    edgeN = perEdgeNormals(V, S);
    vertN = zeros(size(S, 1), 2);
    vertN(S(:, 1), :) = vertN(S(:, 1), :) + edgeN;
    vertN(S(:, 2), :) = vertN(S(:, 2), :) + edgeN;

    N = vertN ./ 2;
    N = normalizerow(N);
end

