
%Calculates straight line going from startPoint to stopPoint with res
%segments
function [V, E] = buildStraightLine(res, startPoint, stopPoint)
    disp =  stopPoint - startPoint ;
    dist = norm(disp);
    dir = disp/dist;
    V = zeros(res+1, 2);
    stepsize = dir/res;
    
    V(:, 1) = linspace(startPoint(1), stopPoint(1), res + 1);
    V(:, 2) = linspace(startPoint(2), stopPoint(2), res + 1);

    E = zeros(res, 2); 
    E(:, 1) = 1:res;
    E(:, 2) = 2:res+1;
    
end
