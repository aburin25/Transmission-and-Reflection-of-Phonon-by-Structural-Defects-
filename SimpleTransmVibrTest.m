function y=SimpleTransmVibrTest(Om, a, A)
% Evaluation of transmission in one dimensional chain with force constant
% a, and force constant defects A. For example the input A=[ 5*a ] suggest
% the presence of a single defect characterized by a force constant 5*a
% within the infinite chain with the force constant a
% Om is the incoming phonon frequency
% ******************************************
%  Output 
% ******************************************
% y.RT(1) - transmission 
% y.RT(2) - reflection
% ******************************************
A=[a a a A a a a]
M=[-Om^2/a+2 -1; 1 0]; 
[V, D]=eig(M); 
D=diag(D); 
disp(D);
T=M; 
for n=1:(max(size(A))-1)
    T=[(-Om^2/A(n+1)+1+A(n)/A(n+1)) -A(n)/A(n+1); 1 0]*T;
end
Res=V^(-1)*T*V;
R=-Res(2,1)/Res(2,2); Tr=Res(1,1)+Res(1,2)*R; 
y.RT=[abs(Tr)^2 abs(R)^2]; 
end