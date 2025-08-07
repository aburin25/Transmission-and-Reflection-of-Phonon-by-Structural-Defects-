function [y, qm, v, z]=SpectraTrL(a, b, r, q)
% the program determines the phonon spectrum of the infinite chain 
% for the fence model and the given wavevector of matching branch
% Vibrations are allowed in two dimensions
% *********************************************************
%  INPUT
% *********************************************************
% a is the nearest neighnor force constant, 
% b is the next neighbor force constant
% r is the tangent of the angle between the nearest neighbor bond 
% and the chain axis   
% Chain period is set to 2
% *********************************************************
%  OUTPUT
% *********************************************************
% y - frequency of acoustic phohon with the wavevector q
% qm - wavevector of acoustic phonon from the other spectral branch 
% possessing the same frequency
% v - group velocity of acoustic phonon
% z - frequency of optical phohon with the wavevector q

a=a/(1+r^2); 
if size(q, 1)==1
    q=q';
end
ompl=a*(1-cos(q))+b*(1-cos(2*q))+r^2*a*(1+cos(q)); 
ommin=a*(1-cos(q))+b*(1-cos(2*q))-r^2*a*(1+cos(q));
V=4*r^2*a^2*sin(q).^2; 
z=sqrt(ompl+sqrt(ommin.^2+V));
y=sqrt(ompl-sqrt(ommin.^2+V));
% group velocity
omplder=a*sin(q)+b*sin(2*q)-a*r^2*sin(q); 
omminder=a*sin(q)+b*sin(2*q)+a*r^2*sin(q);
Vder=4*r^2*a^2*sin(2*q); 
yDer=0.5*(omplder-0.5*(2*ommin.*omminder+Vder)./sqrt(ommin.^2+V))./sqrt(y); 
v=yDer; 
qm=0;  
if size(q, 1)*size(q, 2)==1
    Fact=1; 
    if q<0 
        Fact=-1;
    end
    q=q*Fact; 
    fun1=@(q) a*(1-cos(q))+b*(1-cos(2*q))+r^2*a*(1+cos(q))-sqrt((a*(1-cos(q))+b*(1-cos(2*q))-r^2*a*(1+cos(q))).^2+4*r^2*a^2*sin(q).^2); 
    f0=fun1(q); 
    fun2=@(q) fun1(q)-f0; 
    out1(1) = fzero(fun2,0.001*pi); 
    out1(2)= fzero(fun2,0.9*pi); 
    [~, num]=max(abs(out1-q));
%    disp(out1); 
%    disp(num); 
    qm=acos(cos(out1(num))); 
end
%y=[omac omopt]; 
end