
%Want to know which faces in mesh corresponding to V, F, are interior Faces
%to boundary polygon given by Vboundary, Eboundary
% Inputs:
% Vboundary - vertices for boundary geometry
% Eboundary - edges describing topology of boundary
% V - Vertices of geometry we want to segment
% F - faces of geometry we want to segment
function [intF, extF] = segmentDomain(Vboundary, Eboundary, V, F)
    BC = barycenter(V, F);
    W = winding_number(Vboundary, Eboundary, BC);
    W = (W - min(W))/(max(W) - min(W));
    interiorCondition = W > 0.95;
    extF = F(~interiorCondition, :);
    intF = F(~not(interiorCondition), :);
end