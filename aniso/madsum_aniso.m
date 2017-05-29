function [Eel, Eprm, Eimg, Esum] = madsum_aniso(C0, b0, M, rc, coords, trunc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File: madsum_aniso.m
%%
%% Last Modified : Sun Aug  8 16:09:39 2004
%% Wei Cai, caiwei@stanford.edu
%%
%% References:
%%
%%  1. W. Cai, Atomistic and Mesoscale Modeling of Dislocation Mobility, 
%%     Ph. D. Thesis, Massachusetts Institute of Technology, 2001. 
%%     ( http://asm.mit.edu/caiwei/Download/Thesis/CaiWeiThesis.pdf 
%%       p. 296, for dislocation interaction energy in anisotropic medium. )
%%
%%  2. Wei Cai, Vasily V. Bulatov, Jinpeng Chang, Ju Li and Sidney Yip,
%%     Periodic image effects in dislocation modelling,
%%     Philosophical Magazine A, 83, 539 (2003).
%%     ( for regularization of conditional convergence )
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This program calculates the elastic energy of a dislocation
%  dipole in a supercell under periodic boundary conditions
%   
%  Geometry:
%
%    (x1,y1)   A-------------------------------------
%               \                                    \
%                \      P                             \
%                 \    (-) (x2r,y2r)                   \
%                  \   /                                \
%                   \ /                                  \
%                   (+)-----------------------------------B (x0,0)
%
%  Input:
%     C0:    C(i,j,k,l) general elastic constant tensor (in eV/A^3,  1eV/A^3 = 160.22GPa )
%     b0 = [ bx, by, bz ] Burgers vector in unit cell frame (in Angstrom)
%                        e.g. b0 = [ 1 1 1 ]/2 * 3.1472 (for Molybdenum)
%     M = [ e1, e2, e3 ] coordinate transformation matrix, 
%                        e3 is parallel to dislocation line
%     rc:   cut-off radius
%     coords = [ x0, x1, y1, x2r, y2r ]
%           x0            x coord of B     (0, infty)     (in A)
%           x1            x coord of A     (-infty,infty) (in A)
%           y1            y coord of A     (0,  infty)
%           x2r   reduced x coord of P     (-0.5, 0.5]
%           y2r   reduced y coord of P     (-0.5, 0.5]
%    trunc = [ xcut, ycut ] truncation paramters
%           -xcut:xcut times -ycut:ycut number of images will be added
%     
%  Output:
%     Eel:  total elastic energy, Eel = Eprm + Eimg
%     Eprm: primary dipole energy
%     Eimg: image energy
%     Esum: intermediate results


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step 1: solve the energy prefactor in anisotropic elasticity
%
%Transform C
Mt=M';
for i=1:3, for j=1:3, for k=1:3, for l=1:3,
   tmp=0;
   for ip=1:3, for jp=1:3, for kp=1:3, for lp=1:3,
      tmp=tmp+Mt(i,ip)*Mt(j,jp)*Mt(k,kp)*Mt(l,lp)*C0(ip,jp,kp,lp);
   end; end; end; end;
   C(i,j,k,l)=tmp;
end; end; end; end;

%Transform b
b=(Mt*reshape(b0,3,1))';

%Solve the sextic equation in the coordinate specified by M
[p,A,B,D,poly6,a,c0,c1,c2]=sextic(b,C);

%Calculate interaction matrix h1b(i,n)
% only p and h1b are used in following calculations
% other variables are computed and saved for comparison with earlier codes
% ( ./msumdriver < phRI.dat )
h1=zeros(3,3);
for n=1:3, for i=1:3,
    h1(i,n) = B{n}(i,2,1)*A(1,n)*D(n)+B{n}(i,2,2)*A(2,n)*D(n)...
             +B{n}(i,2,3)*A(3,n)*D(n);
end; end;
h1b=b*h1;
pR=real(p(1:3))';pI=imag(p(1:3))';
hR=real(h1);  hI=imag(h1);
hbR=real(h1b);hbI=imag(h1b);

%Save parameters
phRI=[pR;pI;hR;hI;hbR;hbI;b];
save -ascii -double phRI.dat phRI
disp('interaction matrix saved to file phRI.dat (in eV and A)');
%
%end of Step 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step 2: image summation
%
% (Ref. 1, p. 296)
% W(x,y) = sum_{n=1}^{3} Re[ h1b(n) / (2*pi*I) * ln (x+p_n*y)/rc ]

x0 = coords(1);
x1 = coords(2);
y1 = coords(3);
x2r= coords(4);
y2r= coords(5);

xcut = trunc(1);
ycut = trunc(2);

c1 = [ x0, 0  ]';
c2 = [ x1, y1 ]';
ds = [x2r,y2r ]';
dr = [ c1, c2 ]*ds;

R = norm(dr); dt = dr/R;
prm = log( [dr(1)+dr(2)*p(1:3).'] / rc );
Efac = h1b/(2*pi*sqrt(-1)) *(-0.5);
Eprm = sum(real(Efac.*prm))*(-2);

%Add image contribution (naive summation)
%    image coordinate (x,y)  image strength
image = [ 0,             0,  2
          dr(1),     dr(2), -1
         -dr(1),    -dr(2), -1 ];
       
Esum = [ 0 0 0 ];
for i=-xcut:xcut,
   for j=-ycut:ycut,
      if(i~=0)|(j~=0)
         offset=i*c1+j*c2;
         for k=1:length(image(:,1)),
            dR = offset + image(k,1:2)'; 
            Esum = Esum + image(k,3)*log( [dR(1)+dR(2)*p(1:3).']/rc );
         end
      end
   end
end

%Add ghost contribution (correction for conditional convergence)
%        ghost coordinate (x,y)        ghost strength
ghost = [ c1(1)/2-dr(1)/2,  c1(2)/2-dr(2)/2,  x2r
         -c1(1)/2-dr(1)/2, -c1(2)/2-dr(2)/2, -x2r
          c1(1)/2+dr(1)/2,  c1(2)/2+dr(2)/2, -x2r
         -c1(1)/2+dr(1)/2, -c1(2)/2+dr(2)/2,  x2r 
          c2(1)/2-dr(1)/2,  c2(2)/2-dr(2)/2,  y2r
         -c2(1)/2-dr(1)/2, -c2(2)/2-dr(2)/2, -y2r
          c2(1)/2+dr(1)/2,  c2(2)/2+dr(2)/2, -y2r
         -c2(1)/2+dr(1)/2, -c2(2)/2+dr(2)/2,  y2r ];
         
Egst=[ 0 0 0 ];         
for i=-xcut:xcut,
   for j=-ycut:ycut,
      offset=i*c1+j*c2;
      for k=1:length(ghost(:,1)),
         dR = offset + ghost(k,1:2)';
         Egst = Egst - ghost(k,3)*log( [dR(1)+dR(2)*p(1:3).']/rc );
      end
   end
end

Esum=Esum+Egst;
Eimg = sum(real(Efac.*Esum));
Eel = Eprm + Eimg;
