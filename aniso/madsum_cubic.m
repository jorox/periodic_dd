function [Eel, Eprm, Eimg, Esum] = madsum_cubic(C, b0, M, rc, coords, trunc)
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
%     C  = [C11, C12, C44] cubic elastic constants (in eV/A^3,  1eV/A^3 = 160.22GPa )
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
%

%Prepare Elastic Constants Matrix Cijkl
C11=C(1); C12=C(2); C44=C(3);
C0=zeros(3,3,3,3);
%Cijkl in the compact form
cC0=[
  C11 C12 C12 0   0   0
  C12 C11 C12 0   0   0
  C12 C12 C11 0   0   0
  0   0   0   C44 0   0
  0   0   0   0   C44 0
  0   0   0   0   0   C44
  ];
%Voigt notation
cind=[ 1 6 5 ; 6 2 4 ; 5 4 3 ];
for i=1:3, for j=1:3, for k=1:3, for l=1:3,
   C0(i,j,k,l)=cC0(cind(i,j),cind(k,l));
end; end; end; end;

[Eel, Eprm, Eimg, Esum] = madsum_aniso(C0, b0, M, rc, coords, trunc);
