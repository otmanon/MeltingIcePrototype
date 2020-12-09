function [G,F1,F2] = BIE_line_integral(C,E,N1,N2,eta)

    A = C(E(:,2),:)-C(E(:,1),:);
    % B(i,:,k) = C(E(i,1),:) - V(k,:)
    B = bsxfun(@minus,C(E(:,1),:),permute(eta,[3 2 1]));
    Q = sum(A.*A,2);
    S = sum(B.*B,2);
    R = 2*sum(bsxfun(@times,A,B),2);
    % normal
    %UN = A*[0 -1;1 0];
    %N = normalizerow(UN);
    Anorm = normrow(A);
    BA1 = sum(bsxfun(@times,B,bsxfun(@times,Anorm,N1)),2);

    X = bsxfun(@minus,4*bsxfun(@times,S,Q),R.^2);
    X = max(X,1e-20);
    SRT = sqrt(X);
    %S=S+1e-20;
    L0 = log(S);
    L1 = log(bsxfun(@plus,bsxfun(@plus,S,Q),R));
    A0 = atan2(R,SRT)./SRT;
    A1 = atan2(bsxfun(@plus,2*Q,R),SRT)./SRT;
    A10 = A1-A0;
    L10 = L1-L0;
    
    F1 = bsxfun(@times,BA1/pi,A10);
    F1 = permute(F1,[3 1 2]);
    BC = (C(E(:,2),:)+C(E(:,1),:))/2;
    sqr_d = pdist2(eta,BC).^2;
    [I,J] = find(sqr_d< 2e-7);
    F1(sub2ind(size(F1),I,J)) = -0.5;
    
    if(~isempty(N2))
        BA2 = sum(bsxfun(@times,B,bsxfun(@times,Anorm,N2)),2);
        F2 = bsxfun(@times,BA2/pi,A10);
        F2 = permute(F2,[3 1 2]);
        BC = (C(E(:,2),:)+C(E(:,1),:))/2;
        sqr_d = pdist2(eta,BC).^2;
        [I,J] = find(sqr_d< 2e-7);
        F2(sub2ind(size(F2),I,J)) = -0.5;
    else
        F2=[];
    end
    
    G = ...
    bsxfun(@times,-Anorm/(4*pi), ...
    bsxfun(@times,bsxfun(@minus,4*S,bsxfun(@rdivide,R.^2,Q)),A10)+ ...
    bsxfun(@plus,bsxfun(@times,bsxfun(@rdivide,R,2*Q),L10),L1)-2);
    G = permute(G,[3 1 2]);
    % points lying on edges
    for e = 1:size(E,1) 
        [t,sqr_d] = project_to_lines(eta,C(E(e,1),:),C(E(e,2),:));
        %on_edge = ((abs(sqr_d) < 2e-7) & ((t > -1e-10) & (t < (1+1e-10))));
        on_edge = abs(sqr_d) < 2e-7;
        % edge length
        s = norm(C(E(e,2),:)-C(E(e,1),:));
        % From maple
        toe = s*t(on_edge);
        G(on_edge,e) = ((toe-s).*log((toe-s).^2)-toe.*log(toe.^2)+2*s)/(4*pi);
        on_end = on_edge & (abs(t)<1e-8 | abs(1-t)<1e-8);
        G(on_end,e) = (-s*log(s^2)+2*s)/(4*pi);
    end
   
end

