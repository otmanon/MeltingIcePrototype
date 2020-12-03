%Given a piecewise constant edge motion, find the best possible
% vertex motion that gives rise to that edge motion
%Iputs:
% V - Vertices of geometry
% E - list of edges of geometry
% N - normals defined at each edge of boundary
% M - Midpoint of each edge of boundary
%Outputs:
% VMotion - vertex motion
function VMotion = fitVertexMotion(V, E, N,  EMotion)

    [localV, g2l, l2g, localE] = remove_unreferenced(V, E);

    [numV, dim] = size(localV);
    


    M = fillMNLSE(localV, localE, N);
    A = fillANLSE(localV, localE, N);
    
    M2 = fillMVLSE(localV, localE);
    A2 = fillAVLSE(localV, localE);
    lambda = 0.1;
    I = eye(numV*dim,numV*dim );
    norms = repmat(EMotion,1,size(N,2));
    EdgeMotions =reshape((norms .* N).', 1, [])';
    b = A*EMotion + A2*(EdgeMotions);
    VMotionFlat = (M + lambda*M2)\ b;
    VMotionLocal = vec2mat(VMotionFlat, 2);
    
    VMotion = zeros(size(V, 1), 2);
    VMotion(l2g, :) = VMotionLocal;
end

%%%%%%%%%%%%%%%% Functions to Set up Least Squares %%%%%%%%%%%%%%%%
function M=fillMVLSE(V, E)
    i = [];
    j = [];
    v = [];
    counter = 1;
    %fill M
   [rows cols] = size(E);
    for edgeIndex=1:(rows)
        v1i = E(edgeIndex, 1);
        v2i = E(edgeIndex, 2);
        eLength = edgeLength(V, E, edgeIndex);
        
        i(counter) = 2*v1i - 1;
        j(counter) = 2*v2i - 1;
        v(counter) = eLength/6;
        counter = counter + 1;
        
        i(counter) = 2*v1i;
        j(counter) = 2*v2i;
        v(counter) = eLength/6;
        counter = counter + 1;


        i(counter) = 2* v2i -1;
        j(counter) = 2 * v1i - 1;
        v(counter) = eLength/6;
        counter = counter + 1;
        
        i(counter) = 2*v2i;
        j(counter) = 2*v1i;
        v(counter) = eLength/6;
        counter = counter + 1;

        %%Diagonal elements get 2x larger contribution from this edge
        i(counter) = 2*v1i - 1;
        j(counter) = 2*v1i - 1;
        v(counter) = eLength/3;
        counter = counter + 1;
        
        i(counter) = 2*v1i;
        j(counter) = 2*v1i;
        v(counter) = eLength/3;
        counter = counter + 1;

        i(counter) = 2*v2i - 1;
        j(counter) = 2*v2i - 1;
        v(counter) = eLength/3;
        counter = counter + 1;
        
        i(counter) = 2*v2i;
        j(counter) = 2*v2i;
        v(counter) = eLength/3;
        counter = counter + 1;
    end
    M = sparse(i, j, v);
end

function A=fillAVLSE(V, E)
    i = [];
    j = [];
    v = [];
    %fill A;
    counter = 1;
    [rows cols] = size(E);
    for edgeIndex=1:(rows)
        v1i = E(edgeIndex, 1);
        v2i = E(edgeIndex, 2);
        eLength = edgeLength(V, E, edgeIndex);
        
        i(counter) =2* v1i - 1;
        j(counter) = 2*edgeIndex - 1;
        v(counter) = edgeLength(V, E,edgeIndex)/2;
        counter = counter + 1;
        
        i(counter) = 2*v1i;
        j(counter) = 2*edgeIndex;
        v(counter) = edgeLength(V, E,edgeIndex)/2;
        counter = counter + 1;

        i(counter) = 2*v2i - 1;
        j(counter) = 2*edgeIndex - 1;
        v(counter) = edgeLength(V, E,edgeIndex)/2;
        counter = counter + 1;
        
        i(counter) = 2*v2i;
        j(counter) = 2*edgeIndex;
        v(counter) = edgeLength(V, E,edgeIndex)/2;
        counter = counter + 1;
    end
    A = sparse(i, j, v);
end

function M=fillMNLSE(V, E, Normals)
    i = [];
    j = [];
    v = [];
    counter = 1;
    %fill M
   [rows cols] = size(E);
    for edgeIndex=1:(rows)
        n = Normals(edgeIndex, :);
        nx = n(1);
        ny = n(2);
        nx2 = nx*nx;
        ny2 = ny*ny;
        nxny = nx*ny;
        v1i = E(edgeIndex, 1);
        v2i = E(edgeIndex, 2);
        eLength = edgeLength(V, E, edgeIndex);
        
        %Top Left 4x4 block
        i(counter) = 2*v1i - 1;  
        j(counter) = 2*v1i - 1;
        v(counter) = 2*eLength*nx2/6;
        counter = counter + 1;
        
        i(counter) = 2*v1i - 1;
        j(counter) = 2*v1i;
        v(counter) = 2*eLength*nxny/6;
        counter = counter + 1;
        
        i(counter) = 2*v1i;
        j(counter) = 2*v1i - 1;
        v(counter) = 2*eLength*nxny/6;
        counter = counter + 1;
        
        i(counter) = 2*v1i;
        j(counter) = 2*v1i;
        v(counter) = 2*ny2*eLength/6;
        counter = counter + 1;
        
        %Top right corner        
        i(counter) = 2*v1i - 1;
        j(counter) = 2*v2i - 1;
        v(counter) = eLength*nx2/6;
        counter = counter + 1;
        
        i(counter) = 2*v1i - 1;
        j(counter) = 2*v2i;
        v(counter) = eLength*nxny/6;
        counter = counter + 1;
        
        i(counter) = 2*v1i;
        j(counter) = 2*v2i - 1;
        v(counter) = eLength*nxny/6;
        counter = counter + 1;
        
        i(counter) = 2*v1i;
        j(counter) = 2*v2i;
        v(counter) = eLength*ny2/6;
        counter = counter + 1;
        
        %Bottom left Corner
        i(counter) = 2*v2i - 1;
        j(counter) = 2*v1i - 1;
        v(counter) = eLength*nx2/6;
        counter = counter + 1;
        
        i(counter) = 2*v2i - 1;
        j(counter) = 2*v1i;
        v(counter) = eLength*nxny/6;
        counter = counter + 1;
        
        i(counter) = 2*v2i;
        j(counter) = 2*v1i - 1;
        v(counter) = eLength*nxny/6;
        counter = counter + 1;
        
        i(counter) = 2*v2i;
        j(counter) = 2*v1i;
        v(counter) = eLength*ny2/6;
        counter = counter + 1;
        
        %Bottom Right Corner
        i(counter) = 2*v2i - 1;  
        j(counter) = 2*v2i - 1;
        v(counter) = 2*eLength*nx2/6;
        counter = counter + 1;
        
        i(counter) = 2*v2i - 1;
        j(counter) = 2*v2i;
        v(counter) = 2*eLength*nxny/6;
        counter = counter + 1;
        
        i(counter) = 2*v2i;
        j(counter) = 2*v2i - 1;
        v(counter) = 2*eLength*nxny/6;
        counter = counter + 1;
        
        i(counter) = 2*v2i;
        j(counter) = 2*v2i;
        v(counter) = 2*eLength*ny2/6;
        counter = counter + 1;
        
    end
    M = sparse(i, j, v);
end

function A=fillANLSE(V, E, Normals)
    i = [];
    j = [];
    v = [];
    %fill A;
    counter = 1;
    [rows cols] = size(E);
    for edgeIndex=1:(rows)
        n = Normals(edgeIndex, :);
        nx = n(1);
        ny = n(2);
        v1i = E(edgeIndex, 1);
        v2i = E(edgeIndex, 2);
        eLength = edgeLength(V, E, edgeIndex);
        
        i(counter) = 2* v1i - 1;
        j(counter) = edgeIndex;
        v(counter) = eLength*nx/2;
        counter = counter + 1;
        
        i(counter) = 2*v1i;
        j(counter) = edgeIndex;
        v(counter) = eLength*ny/2;
        counter = counter + 1;

        i(counter) = 2*v2i - 1;
        j(counter) = edgeIndex ;
        v(counter) = eLength*nx/2;
        counter = counter + 1;
        
        i(counter) = 2*v2i;
        j(counter) =  edgeIndex;
        v(counter) = eLength*ny/2;
        counter = counter + 1;
    end
    A = sparse(i, j, v);
end

%%%%%%%%%%%%%%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%

%Gets edge length of given edge
function e=edgeLength(V, E, index) 
    v1i = E(index, 1);
    v2i = E(index, 2);
    
    v1 = V(v1i, :);
    v2 = V(v2i,:);
    
    e = norm(v2 - v1);
end



function mat=vec2mat(vec, dim)
    mat = zeros(length(vec)/2, dim);
    for index=1:length(vec)/2
        mat(index, :) = [vec(2*index-1) vec(2*index)];
    end
end
function drawLines(V1, V2, color)
    [rows cols] = size(V1);
    for index=1:rows
        line([V1(index, 1) V2(index, 1)], [V1(index, 2) V2(index, 2)], 'Color',color, 'LineWidth',2)
    end
end
