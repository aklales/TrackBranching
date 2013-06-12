clear all;



DontReload = 0;

Propagate_interpolate

density = betterbinz(1024,256,xx,yy,Lx,Ly);

DontReload = 1;

for ii = 1:10
    ii
    Propagate_interpolate
    currentbins = betterbinz(1024,256,xx,yy,Lx,Ly);
    density = density + currentbins;
end