%test madsum_iso

mu=1; nu=0.3; b=[1 2 3]; rc=1; coord=[1 0.2 2 0 0.5]; cut=[10 10];

[Eel, Eprm, Eimg, Esum] = madsum_iso(mu, nu, b, rc, coord, cut);

Eel,Eprm,Eimg,Esum'

