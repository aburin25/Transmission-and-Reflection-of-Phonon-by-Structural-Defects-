function v=GroupVel(a, b, r, q)
% Evaluates group velocities for the given wavevector q within the fence
% model 
ompl=a*(1-cos(q))+b*(1-cos(2*q))+r^2*a*(1+cos(q)); 
ommin=a*(1-cos(q))+b*(1-cos(2*q))-r^2*a*(1+cos(q));
V=4*r^2*a^2*sin(q).^2; 
z=sqrt(ompl+sqrt(ommin.^2+V));
y=sqrt(ompl-sqrt(ommin.^2+V));
% group velocity
omplder=a*sin(q)+2*b*sin(2*q)-a*r^2*sin(q); 
omminder=a*sin(q)+2*b*sin(2*q)+a*r^2*sin(q);
Vder=4*r^2*a^2*sin(2*q); 
yDer=0.5*(omplder-0.5*(2*ommin.*omminder+Vder)./sqrt(ommin.^2+V))./y; 
v=yDer;
end