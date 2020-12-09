function N = perEdgeNormals(C, E);    
    A = C(E(:,2),:)-C(E(:,1),:);
    UN = A*[0 -1;1 0];
    N = normalizerow(UN);
end