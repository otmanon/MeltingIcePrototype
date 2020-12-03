% Calculates flux along boundary given by S.
% Input:
% V - Vertices of Geometry
% F - Faces of geometry
% Outputs:
% dTdn - flux in normal direction of boundary
% M     - midpoint of boundary edges
function [dTdn, N, M, S] = getFluxAlongBoundary(V, intF, F, TGrad)
    S = boundary_faces( intF);
    M = midpoints(V, S);
    N = TwoDNormals(V, S);
    ET = edge_triangle_adjacency(F, S); %edge triangle adjacency
    
    %Get Barycenter of each face
    %Need the index of interior face and exterior face. 
    %Ensure each face in interior is really an interior face by checking
    %with normal
    %Take gradient entries indexed by interior faces, subtract from them
    %gradient entries indexed by exterior faces. Take resulting vector and
    %dot with normal.
    BC = barycenter(V, F);
    FS = ET(:, 1); %solid triangle
    FL = ET(:, 2); % liquid triangle
    for i = 1:size(S, 1)
        if(FL(i) > 0)
            lPos = BC(FL(i), :);
            mid = M(i, :);
            n = N(i, :);
            sign = dot(n, lPos - mid);
            if(sign < 0)
                %swap FS with FL
                temp = FL(i);
                FL(i) = FS(i);
                FS(i) = temp;
            end
        end
       
    end
    
    TGradS = TGrad(FS, :);
    FLValidInd = find(FL > 0); %Some values of Liquid boundary will be -1.
    TGradL = zeros(size(TGradS, 1), 2);
    TGradL(FLValidInd, :) = TGrad(FL(FLValidInd), :);
    
    dTdn = TGradS - TGradL;
    dTdn = dot(dTdn, N, 2);
    
    
end

function M = midpoints(V, S);
    V1 = V(S(:, 1), :);
    V2 = V(S(:, 2), :);
    M = V1 + V2;
    M = M ./ 2;
end

function N = TwoDNormals(V, S)
    V1 = V(S(:, 1), :);
    V2 = V(S(:, 2), :);
    D(:, 1) = (V2(:, 1) - V1(:, 1));
    D(:, 2) = (V2(:, 2) - V1(:, 2));
    Dnorm = vecnorm(D,2, 2);
    D = D ./ Dnorm;
    N = zeros(length(D), 2);
    N(:, 1) = D(:, 2);
    N(:, 2) = -D(:, 1);
end
