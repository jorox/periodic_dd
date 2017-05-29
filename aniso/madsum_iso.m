function [Eel, Eprm, Eimg, Esum] = madsum_iso(mu, nu, b, rc, coords, trunc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File: madsum_iso.m
%%
%% Last Modified : Sun Aug  8 16:09:39 2004
%% Wei Cai, caiwei@stanford.edu
%%
%% References:
%%
%%     Wei Cai, Vasily V. Bulatov, Jinpeng Chang, Ju Li and Sidney Yip,
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
%     mu:   shear modulus (in eV/A^3,  1eV/A^3 = 160.22GPa )
%     nu:   Poisson's ratio
%     b      = [ bx, by, bz ]              Burgers vector
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
%     Esum = [ lnR, xxR2, yyR2, xyR2 ] : intermediate results
%            lnR :    sum of ln(R/rc)
%            xxR2:    sum of x*x/R^2
%            yyR2:    sum of y*y/R^2
%            xyR2:    sum of x*y/R^2
%

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
prm = [log(R/rc), dt(1)*dt(1), dt(2)*dt(2), dt(1)*dt(2)];
Efac = -(mu/4/pi)*[ (b(1)*b(1)+b(2)*b(2))/(1-nu)+b(3)*b(3), ...
        b(2)*b(2)/(1-nu), b(1)*b(1)/(1-nu), -2*b(1)*b(2)/(1-nu) ];
Eprm = sum(Efac.*prm)*(-2);

%Add image contribution (naive summation)
%    image coordinate (x,y)  image strength
image = [ 0,             0,  2
          dr(1),     dr(2), -1
         -dr(1),    -dr(2), -1 ];
       
Esum = [ 0 0 0 0 ];
for i=-xcut:xcut,
   for j=-ycut:ycut,
      if(i~=0)|(j~=0)
         offset=i*c1+j*c2;
         for k=1:length(image(:,1)),
            dR = offset + image(k,1:2)'; R = norm(dR); dt = dR/R;
            Esum = Esum + image(k,3)*[log(R/rc), dt(1)*dt(1), dt(2)*dt(2), dt(1)*dt(2)];
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
         
Egst=[ 0 0 0 0 ];         
for i=-xcut:xcut,
   for j=-ycut:ycut,
      offset=i*c1+j*c2;
      for k=1:length(ghost(:,1)),
         dR = offset + ghost(k,1:2)'; R = norm(dR); dt = dR/R;
         Egst = Egst - ghost(k,3)*[log(R/rc), dt(1)*dt(1), dt(2)*dt(2), dt(1)*dt(2)];
      end
   end
end

Esum=Esum+Egst;
Eimg = sum(Efac.*Esum);
Eel = Eprm + Eimg;
