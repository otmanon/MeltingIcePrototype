function [W, psi] = solveBEM(C, S, intS, extS, Tb1, Tb2, V)
    N = perEdgeNormals(C, S);
   
    %Tb1 and Tb2 are vertex based... gotta make theme edge based

    
    
    BC = barycenter(C, S);
    BCint = barycenter(C,intS);
    bint = unique(intS);
    bext = unique(extS);
    N1int = N(bint, :); N2int = -N(bint, :);
    N1ext = N(bext, :); N2ext = -N(bext, :);
    
    Teb = Tb1(intS(:, 1)) + Tb1(intS(:, 2));
    Teb = Teb ./ 2;
    Teb = [Teb; Tb1(end -3:end)];
    Tb1 = Teb;
    Tb2 = Teb;
    
    [Gint,Fint1,Fint2]=BIE_line_integral(C,intS,N1int,N2int,BCint);
    [Gvint,Fv1int,Fv2int]=BIE_line_integral(C,intS,N1int,N2int,V);        
    psiInt = Gint\(Fint1*Tb1(1:end-4)+Fint2*Tb2(1:end-4));
    
    [Gext,Fext1,Fext2]=BIE_line_integral(C,S,N,-N,BC);
    [Gvext,Fv1ext,Fv2ext]=BIE_line_integral(C,S,N, -N,V);
    psiExt = Gext\(Fext1*Tb1+Fext2*Tb2);
    
    
    psi =  psiInt- psiExt(bint); 
    
    W = -Gvint*psiInt+Fv1int*Tb1(bint)+Fv2int*Tb2(bint);

    W = W  - Gvext*psiExt+Fv1ext*Tb1+Fv2ext*Tb2;



end

