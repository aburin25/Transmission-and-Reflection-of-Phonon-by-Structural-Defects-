function y=TransfMatrGen(A, B, r, n, Om, mm)
b=B; 
% A - column vector of nearest neighbor force constants of the length N+2,
% A(n) is the force constant between sites n-1 and n
% B - column vector of next neighbor force constants of the length N+2, 
% B(n) is the force constant between sites n-2 and n
% m - column vector of masses of the length N+1,
% r - column vector of tangents of angles between the molecular axis and the n-1 - n bond of the length N+1,
% n -integer number from 1 to N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outcome y - transfer matrix 7x7 shifting  column vector [ux(n+1), ux(n),
% ux(n-1), ux(n-2), uy(n+1), uy(n), uy(n-1)] by one towards increasing n
disp(n); 
disp(size(A)); disp(size(b)); disp(size(r)); disp(size(mm)); 
M1=[-A(n+1)/b(n+2) -mm(n)*Om^2/(b(n+2))+(A(n+1)+A(n))/b(n+2)+1+b(n)/b(n+2) -A(n)/b(n+2) -b(n)/b(n+2) A(n+1)*r(n+1)/b(n+2) 1/b(n+2)*(r(n+1)*A(n+1)-r(n)*A(n)) -r(n)*A(n)/b(n+2); eye(3) zeros(3, 4); zeros(3, 4) eye(3)];
M2=[eye(4) zeros(4, 3); 1/r(n+2) -1/r(n+2)+r(n+1)*A(n+1)/(A(n+2)*r(n+2)^2) -A(n+1)*r(n+1)/(r(n+2)^2*A(n+2)) 0 mm(n+1)*Om^2/(A(n+2)*r(n+2)^2)-(1+A(n+1)*r(n+1)^2/(A(n+2)*r(n+2)^2)) -A(n+1)*r(n+1)^2/(A(n+2)*r(n+2)^2) 0; zeros(2, 4) eye(2) zeros(2, 1)];
y=M2*M1;
end

% M1=[-a/b -Om^2/b+2*a/b+2  -a/b -1  a*r/b 0 -a*r/b; eye(3) zeros(3, 4); zeros(3, 4) eye(3)];
%  M2=[eye(4) zeros(4, 3); 1/r 0 -1/r 0 Om^2/(a*r^2)-2 -1 0; zeros(2, 4) eye(2) zeros(2, 1)]; 
 