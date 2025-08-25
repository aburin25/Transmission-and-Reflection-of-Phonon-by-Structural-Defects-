function y=TransmChain(a, b, r, m, A, B, R, MM, Om)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluation of transmission and reflection of longitudinal wave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a is the nearest neighbor force constant in ordered chain (Sample a=4)
% b is the next neighbor force constant  in ordered chain (b=1)
% r is the tangent of angle between the molecular axis (x) and bond
% direction in ordered chain (r=1)
% m is the mass of the unit cell (m=1)
% A - column vector of nearest neighbor force constants of the length N+2
% for chain including defects (A=[4 4]')
% A(n) is the force constant between sites n-1 and n
% B - column vector of next neighbor force constants of the length N+2, 
% B(n) is the force constant between sites n-2 and n (B=[2 1]')
% M - column vector of masses of the length N+1 (M=[1 1])
% R - column vector of tangents of angles between the molecular axis and the n-1 - n bond of the length N+1,
% domain (R=[1 1])
% Input corresponds to a single defect of the nearest neighbor force
% constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y.RT - 4x1 column vector 
% y.RT(1) - reflection of longitudinal to longitudinal
% y.RT(2) - reflection of longitudinal to transverse
% y.RT(3) - trabnsmission of longitudinal to longitudinal
% y.RT(4) - trabnsmission of longitudinal to transverse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=a/(1+r^2); A=A./(1+R.^2); 

 M1=[-a/b -m*Om^2/b+2*a/b+2  -a/b -1  a*r/b 0 -a*r/b; eye(3) zeros(3, 4); zeros(3, 4) eye(3)];
 M2=[eye(4) zeros(4, 3); 1/r 0 -1/r 0 m*Om^2/(a*r^2)-2 -1 0; zeros(2, 4) eye(2) zeros(2, 1)]; 
 M=M2*M1; %Msp=M2sp*M1sp; 
 [V, D]=eig(M); 
 % exclude zero 
 D=diag(D); 
 [~, k]=min(abs(D));
 if k<7
    D=D([1:(k-1) (k+1):7 k]); 
    V=V(:, [1:(k-1) (k+1):7 k]); 
 end
 [~, k]=min(abs(D(1:6)));
 if k<6
    D=D([1:(k-1) (k+1):6 k 7]); 
    V=V(:, [1:(k-1) (k+1):6 k 7]); 
 end
  [~, k]=min(imag(D(1:5)));
 if k<5
    D=D([1:(k-1) (k+1):5 k 6:7]); 
    V=V(:, [1:(k-1) (k+1):5 k 6:7]); 
 end
 [~, k]=max(imag(D(1:4)));
 if k<4
    D=D([1:(k-1) (k+1):4 k 5:7]); 
    V=V(:, [1:(k-1) (k+1):4 k 5:7]);
 end
 q1=pi-asin(imag(D(4)));
 %disp(q1); 
 v1=abs(GroupVel(a, b, r, q1)); 
 %disp(v1); 
 [~, k]=min(imag(D(1:3)));
 if k<3
    D=D([1:(k-1) (k+1):3 k 4:7]); 
    V=V(:,[1:(k-1) (k+1):3 k 4:7]); 
 end
 [~, k]=max(imag(D(1:2)));
 if k<2
    D=D([2 1 3:7]); 
    V=V(:, [2 1 3:7]); 
 end
 q=asin(imag(D(2)));
 v=GroupVel(a, b, r, q);
 %disp(v);
 for ii=2:5
    V(:, ii)=V(:, ii)/norm(V(:, ii)); 
    V(:, ii)=V(:, ii)/sqrt(abs(V(1, ii))^2+abs(V(5, ii))^2); 
 end
 %y.Bas=V; 
 %VV=V; 
 %D(k)=[]; VV(:, k)=[];
 %Ans=zeros(7, 6);
 AA=[a a A' a a];
%ode rr=[r r r r r r r*0.5 r*2 r r r r r]; 
  rr=[r r R' r r]; 
 B=[b b B' b b]; 
 MM=[m m MM' m m];
 T=M; 
 for ii=1:(max(size(AA))-2)
     disp(ii); 
%    T1=TransfMatrrrr(AA, bb, rr, ii, Om,MM);
    T1=TransfMatrGen(AA, B, rr, ii, Om, MM); 
    disp(norm(T1-M)); 
    T=T1*T; 
 end
 % Find and expand transmitted vector over those eigenvectors
 for k=1:7
     disp(D(k));
%    Vin=VV(:, k)*D(k)^(-5); 
%    Vfin=zeros(size(V, 1), 4); 
%    Vfin1=M^4*Msp*M^5*Vin; 
    % Vfin2=M*Vfin1; 
    % Vfin3=M*Vfin2;
    % Vfin4=M*Vfin3; 
    % Vfin5=M*Vfin4;
    % Vfin6=M*Vfin5; 
    % Vfin=[Vfin1 Vfin2 Vfin3 Vfin4 Vfin5 Vfin6]; 
%    Ans(:, k)=T*Vin;
  
 end
   Ans=T*V; 
 T=V^(-1)*Ans; 
 % Ans(7, :)=[]; 
 % DD=repmat(D, [1, 6]); 
 % T=Ans.*DD.^(-5);
 %y.TrM=T;
  T(7, :)=[]; 
 VVV=-T(1:5, 2);
 VVV1=-T(1:6, 2);
 MMM2=[T(1, 1) T(1, 3) T(1, 4) 0 0; T(2, 1) T(2, 3) T(2, 4) -1 0; T(3, 1) T(3, 3) T(3, 4) 0 0; T(4, 1) T(4, 3) T(4, 4) 0 0; T(5, 1) T(5, 3) T(5, 4) 0 -1];
% MMM=[T(:, [1 3 4]) [0 0 0;   -1 0 0; 0 0 0; 0 0 0; 0 -1 0; 0 0 -1]];
% MMM1=[T(:, [1 3 5]) [0 0 0; -1 0 0; 0 0 0; 0 -1 0; 0 0 0; 0 0 -1]];
 RT=abs(MMM2^(-1)*VVV).^2; 
 RT(3)=RT(3)*v1/v;
 RT(5)=RT(5)*v1/v;
 y.RT=RT(2:5); 
 %y.RT1=abs(MMM1^(-1)*VVV1).^2; 
 % y.Tr=y.Rt(4); 
 % y.Ref1=y.RT(2);
 % y.Ref2=y.RT(3);
end