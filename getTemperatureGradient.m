%Takes mesh and temperatures stored at vertices and returns temperatuire gradient
function TGrad = getTemperatureGradient(V, F, T)
    G = grad(V, F);
    TGrad = reshape(G*T,size(F,1),size(V,2));
end

