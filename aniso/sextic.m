function [p,A,B,D,Poly6,a,c0,c1,c2]=sextic(b,C)
%classical sextic formalism of anisotropic linear elasticity for
%dislocations
%
%[p,A,B,D]=sextic(b,C)
%Inputs:  b(3)          Burger's vector
%         C(3,3,3,3)    Elastic Constant Matrix
%Outputs: p(6)          Six roots of the sextic polynomial p(n)
%         A(3,3)        A(k,n)
%         B{3}(3,3,3)   B{n}(i,j,k)
%         D(3)          D(n)
%
%Displacement field of dislocation (along z direction)
%         u(k)=Re(-1/(2*pi*i)*Sum_{n=1}^{3}(A(k,n)*D(n)*ln(eta(n)))) ,
%         where eta(n)=x+p(n)*y

%Construct Poly6
c0=zeros(3,3);c1=c0;c2=c0;
c0(:,:)=C(:,1,:,1);
c1(:,:)=C(:,1,:,2)+C(:,2,:,1);
c2(:,:)=C(:,2,:,2);

Poly6=zeros(7,1);
Poly6(7)=det(c0);
Poly6(6)...
    =det([c0(:,1),c0(:,2),c1(:,3)])...
    +det([c0(:,1),c1(:,2),c0(:,3)])...
    +det([c1(:,1),c0(:,2),c0(:,3)]);
Poly6(5)...
    =det([c0(:,1),c0(:,2),c2(:,3)])...
    +det([c0(:,1),c2(:,2),c0(:,3)])...
    +det([c2(:,1),c0(:,2),c0(:,3)])...
    +det([c0(:,1),c1(:,2),c1(:,3)])...
    +det([c1(:,1),c0(:,2),c1(:,3)])...
    +det([c1(:,1),c1(:,2),c0(:,3)]);
Poly6(4)...
    =det([c0(:,1),c1(:,2),c2(:,3)])...
    +det([c1(:,1),c0(:,2),c2(:,3)])...
    +det([c0(:,1),c2(:,2),c1(:,3)])...
    +det([c1(:,1),c2(:,2),c0(:,3)])...
    +det([c2(:,1),c0(:,2),c1(:,3)])...
    +det([c2(:,1),c1(:,2),c0(:,3)])...
    +det([c1(:,1),c1(:,2),c1(:,3)]);
Poly6(3)...
    =det([c1(:,1),c1(:,2),c2(:,3)])...
    +det([c1(:,1),c2(:,2),c1(:,3)])...
    +det([c2(:,1),c1(:,2),c1(:,3)])...
    +det([c0(:,1),c2(:,2),c2(:,3)])...
    +det([c2(:,1),c0(:,2),c2(:,3)])...
    +det([c2(:,1),c2(:,2),c0(:,3)]);
Poly6(2)...
    =det([c1(:,1),c2(:,2),c2(:,3)])...
    +det([c2(:,1),c1(:,2),c2(:,3)])...
    +det([c2(:,1),c2(:,2),c1(:,3)]);
Poly6(1)=det(c2);

%Solve the Polynomial
cp=roots(Poly6);
%cp,Poly6

%Sort roots
[tmp,ind]=sort(-imag(cp));
p=cp;
p(1:3)=cp(ind(1:3));
p(4:6)=conj(p(1:3));

%p
%pause
%Construct Matrix A
a=cell(3,1);
I=sqrt(-1);
if abs(p(1)-I)<1e-4 & abs(p(2)-I)<1e-4 & abs(p(3)-I)<1e-4
  disp('degeneracy!');
  for n=1:3,
    a{n}=c0+c1*p(n)+c2*p(n)^2;
  end
  A=[1 -I 0; I 1 0; 0 0 1];
else
  A=zeros(3,3);
  for n=1:3,
    a{n}=c0+c1*p(n)+c2*p(n)^2;
    [v,d]=eig(a{n});
    [s,t]=sort(abs(diag(d)));
    A(:,n)=v(:,t(1));
  end
end
%a{1},a{2},a{3},A,a{1}*A,a{2}*A,a{3}*A
%pause
%Construct Matrix B
B=cell(3,1);
for n=1:3,
  B{n}=C(:,:,:,1)+C(:,:,:,2)*p(n);
end

F=zeros(3,3);
for j=1:3,
  for n=1:3,
    F(j,n)=B{n}(j,2,1)*A(1,n)+B{n}(j,2,2)*A(2,n)+B{n}(j,2,3)*A(3,n);
  end
end

G=zeros(6,6);rhs=zeros(6,1);
G(1:3,1:3)=real(A);
G(1:3,4:6)=-imag(A);
G(4:6,1:3)=real(F);
G(4:6,4:6)=-imag(F);
rhs(1:3)=b;
lhs=G\rhs;

%A,F,G,rhs,lhs
%pause
D=lhs(1:3)+sqrt(-1)*lhs(4:6);
