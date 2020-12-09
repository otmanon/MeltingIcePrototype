function [W, psi] = solveBEM(C, S, intS, extS, Tb1, Tb2, V)
    N = perEdgeNormals(C, S);
   
    BC = barycenter(C, S);
    BCint = barycenter(C,intS);
    bint = unique(intS);
    bext = unique(extS);
    N1int = N(bint, :); N2int = -N(bint, :);
    N1ext = N(bext, :); N2ext = -N(bext, :);
    
    
    [Gint,Fint1,Fint2]=BIE_line_integral(C,intS,N1int,N2int,BCint);
    [Gvint,Fv1int,Fv2int]=BIE_line_integral(C,intS,N1int,N2int,V);        
    psiInt = Gint\(Fint1*Tb1(bint)+Fint2*Tb2(bint));
    
    [Gext,Fext1,Fext2]=BIE_line_integral(C,S,N,-N,BC);
    [Gvext,Fv1ext,Fv2ext]=BIE_line_integral(C,S,N, -N,V);
    psiExt = Gext\(Fext1*Tb1+Fext2*Tb2);
    
    
    psi = psiInt - psiExt(bint);
    
    W = -Gvint*psiInt+Fv1int*Tb1(bint)+Fv2int*Tb2(bint);

    W =  - Gvext*psiExt+Fv1ext*Tb1+Fv2ext*Tb2;



end

