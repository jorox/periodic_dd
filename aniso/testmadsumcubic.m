%test madsum_cubic

%b0= [1 1 1]/2 * 3.1472;                  %edge dislocation for Mo in Angstrom
%M = [ 1/sqrt(6)   1/sqrt(2)   1/sqrt(3)  %for Mo <111> screw
%      1/sqrt(6)  -1/sqrt(2)   1/sqrt(3)
%     -2/sqrt(6)   0           1/sqrt(3) ];
%M = [ 1/sqrt(3)  -1/sqrt(2)  1/sqrt(6)   %for Mo <110> edge b=<111>
%      1/sqrt(3)   0         -2/sqrt(6) 
%      1/sqrt(3)   1/sqrt(2)  1/sqrt(6) ];

b0=[-1 1 0]/2 * 5.4310;                     %screw for Si in Angstrom
M = [ 1/sqrt(3)   1/sqrt(6)  -1/sqrt(2)     %for <111> shuffle screw b=<110>/2
      1/sqrt(3)   1/sqrt(6)   1/sqrt(2)   
      1/sqrt(3)  -2/sqrt(6)   0          ];

%Elastic Constants of Mo (FS-EAM potential, documented)
%C11=4.647e11; C12=1.615e11; C44=1.089e11; %in Pa (99.1% for screw)

%Elastic Constants of Si (SW potential, documented)
C11=1.616e11; C12=0.816e11; C44=0.603e11; %in Pa (97.07% agreement)


C=[C11,C12,C44]/1.6022e11; %(convert into eV/A^3)

rc=3.84; coord=[3*sqrt(3)*5.431 0 8mill*sqrt(6)*5.431 0 0.5]; cut=[15 15];

[Eel, Eprm, Eimg, Esum] = madsum_cubic(C, b0, M, rc, coord, cut);
Eel,Eprm,Eimg,Esum'
