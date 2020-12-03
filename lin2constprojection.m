clc; 
numNodes = 3;
numEdges = 3;
visual_scale = 0.2;
%V = [0 0; 1 0; 1 1]%; 3 0; 4 0; 5 0; 6 0; 7 0; 8 0; 9 0; 10 0];
%E = [ 1 2; 2 3; 3 1]%;  3 4; 4 5; 5 6; 6 7; 7 8; 8 9; 9 10; 10 11];

V = [0 0; 0.5 0.05; 1 0.05; 1 0.5; 1 1; 0.5 1; 0 1; 0 0.5];% ; 3 0; 4 0; 5 0; 6 0; 7 0; 8 0; 9 0; 10 0];
E = [ 1 2;  2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 1] ;% 3 4; 4 5; 5 6; 6 7; 7 8; 8 9; 9 10; 10 11;];
V = V
[numV, dim] = size(V);
[rows cols] = size(E);
midP = zeros(rows, 2);
normals = zeros(rows, 2);
counter = 1;

%Calculates normals, and midpoints for visualization
for index=1:rows
    v1i = E(index, 1);
    v2i = E(index, 2);
    
    v1 = V(v1i, :);
    v2 = V(v2i,:);
    
    disp = v2 - v1;
    midP(counter, :) = [(v1(1) + v2(1))/2 (v1(2) + v2(2))/2];
    
    normal =  [disp(2) -disp(1)];
    normals(counter, :) = -(normal)/norm(normal);
    
    counter = counter + 1;

end

edgesVp = normals;
%edgesVp(:, 2) = 3*edgesVp(:, 2) + 3*(rand(rows, 1) - 1)/2
edgesVpn = projectToNormal(edgesVp, normals);
edgesVp = transpose(reshape(edgesVp.', 1, []));

%[C, M] = fd_laplacian(V, E);
[C, M] = fd_operators(V);
N = perVertexNormals(V, E, normals);
LV = inv(M)*C*V;

curvature = transpose(dot(transpose(LV), transpose(N))); %project onto per vertex normals
curvatureNormal = LV(:, :)./curvature;


drawMesh(V, E);

drawVertices(V);

%M = fillMVLSE(V, E)%, normals)
%A = fillAVLSE(V, E)%, normals)

M = fillMNLSE(V, E, normals)
A = fillANLSE(V, E, normals)

lambda = 0.01;
I = eye(numV*dim,numV*dim );
b = A*edgesVpn;
edgesV = (M - lambda*I)\ b
edgesV = vec2mat(edgesV, 2)
edgesVp = vec2mat(edgesVp, 2)

interpP = interpolateVector(V, 30, V, E);
interpV = interpolateVector(edgesV, 30, V, E);



drawLines(interpP, interpP + interpV*visual_scale, 'green')
drawLines(midP, midP + edgesVp*visual_scale, 'red')
%drawLines(V, V + edgesV*visual_scale, 'blue')

drawLines(V, V + N*visual_scale, 'blue')
drawLines(V, V + perVertexNormal*visual_scale, 'yellow')
%%%%%%%%%%%%%%%% Functions to Set up Least Squares %%%%%%%%%%%%%%%%
function M=fillMVLSE(V, E)
    i = [];
    j = [];
    v = [];
    counter = 1;
    %fill M
   [rows cols] = size(E)
    for edgeIndex=1:(rows)
        v1i = E(edgeIndex, 1);
        v2i = E(edgeIndex, 2);
        eLength = edgeLength(V, E, edgeIndex)
        
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
    [rows cols] = size(E)
    for edgeIndex=1:(rows)
        v1i = E(edgeIndex, 1);
        v2i = E(edgeIndex, 2);
        eLength = edgeLength(V, E, edgeIndex)
        
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
   [rows cols] = size(E)
    for edgeIndex=1:(rows)
        n = Normals(edgeIndex, :);
        nx = n(1);
        ny = n(2);
        nx2 = nx*nx;
        ny2 = ny*ny;
        nxny = nx*ny;
        v1i = E(edgeIndex, 1);
        v2i = E(edgeIndex, 2);
        eLength = edgeLength(V, E, edgeIndex)
        
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
    [rows cols] = size(E)
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
function N=perVertexNormals(V, E, edgeNormals)

    [rows cols] = size(V);
    N = zeros(rows, cols);
    W = zeros(rows, 1)
    for index = 1:rows
        v1i = E(index, 1);
        v2i = E(index, 2);
        l = edgeLength(V, E, index);
        N(v1i, : ) = N(v1i, : ) + edgeNormals(index, :)*l;
        N(v2i, : ) = N(v2i, : ) + edgeNormals(index, :)*l;
        W(v1i) = W(v1i) + l;
        W(v2i) = W(v2i) + l;
    end
    N = N./W;
    N = N ./ transpose(vecnorm(transpose(N)));
end
%Gvien a list of vectors, and a list of normals of the same size, return a
%list of the vector projected on the normal
function out=projectToNormal(vectors, normals)
out = [];
[rows cols] = size(vectors);
    for index=1:rows
        out(index, 1) = dot(vectors(index, :), normals(index, :)) ;
    end
end
function r=interpolateVector(vec, numSamples, V, E)
    [rows, cols] = size(vec);
    r = zeros(rows*(numSamples-1), cols);
    counter = 1;
    [rows cols] = size(E);
    for edgeIndex=1:(rows)
        v1i = E(edgeIndex, 1);
        v2i = E(edgeIndex, 2);
        v1 = vec(v1i, :);
        v2 = vec(v2i, :);
        
        eLength = edgeLength(V, E, edgeIndex);
        s = 0;
        for sampleNum = 1:numSamples-1
            s = s + 1/numSamples;
            mid = s * v2 + (1 - s) * v1;
            r(counter, :) = mid;
            counter = counter + 1;
        end
        
    end
end

%Gets edge length of given edge
function e=edgeLength(V, E, index) 
    v1i = E(index, 1);
    v2i = E(index, 2);
    
    v1 = V(v1i, :);
    v2 = V(v2i,:);
    
    e = norm(v2 - v1);
end

function drawMesh(V, E)
    [rows cols] = size(E);
    for index=1:(rows)
        v1i = E(index, 1);
        v2i = E(index, 2);

        v1 = [V(v1i, 1) V(v1i, 2)];
        v2 = [V(v2i, 1) V(v2i, 2)];
        axis equal
        line([v1(1)  v2(1)], [v1(2)  v2(2)], 'Color','black', 'LineWidth',3)
    end
end

function drawVertices(V)
    hold on
    plot(V(:, 1), V(:, 2), 'r*')
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
