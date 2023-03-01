function ThetaInDegrees = vectorAngles(u,v)

% https://www.mathworks.com/matlabcentral/answers/501449-angle-betwen-two-3d-vectors-in-the-range-0-360-degree

% the order of u and v matters: u represent the reference vec (positive x axis), so the
% positive/negative sign of theta depends on the angle relationship of v toward u

% CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
% ThetaInDegrees = real(acosd(CosTheta));
n=[0 0 1];

if length(u)==2
    u=[u,0];
end
if length(v)==2
    v=[v,0];
end

x = cross(u,v);
c = sign(dot(x,n)) * norm(x);
ThetaInDegrees = atan2d(c,dot(u,v));