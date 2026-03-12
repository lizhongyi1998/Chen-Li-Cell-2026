function [fpx,fpy,fpz,varb]=func_p(vx,vy,vz,px,py,pz,const,varb)

dx=const.dx;
dy=const.dy;
dz=const.dz;
epsl=const.epsl;
gamap=const.gamap;

dxpx=varb.dxpx;
dxpy=varb.dxpy;
dxpz=varb.dxpz;
dypx=varb.dypx;
dypy=varb.dypy;
dypz=varb.dypz;
dzpx=varb.dzpx;
dzpy=varb.dzpy;
dzpz=varb.dzpz;

hx=varb.hx;
hy=varb.hy;
hz=varb.hz;

dyvx=DiffY3D(vx,dy);
dxvy=DiffX3D(vy,dx);
dzvy=DiffZ3D(vy,dz);
dyvz=DiffY3D(vz,dy);
dzvx=DiffZ3D(vx,dz);
dxvz=DiffX3D(vz,dx);

Dxx=DiffX3D(vx,dx);
Dyy=DiffY3D(vy,dy);
Dzz=DiffZ3D(vz,dz);
Dxy=(dyvx+dxvy)/2;
Dyz=(dzvy+dyvz)/2;
Dxz=(dzvx+dxvz)/2;
Oxy=(dyvx-dxvy)/2;
Oyz=(dzvy-dyvz)/2;
Oxz=(dzvx-dxvz)/2;

varb.Dxx=Dxx;
varb.Dyy=Dyy;
varb.Dzz=Dzz;

Dpx=Dxx.*px+Dxy.*py+Dxz.*pz;
Dpy=Dxy.*px+Dyy.*py+Dyz.*pz;
Dpz=Dxz.*px+Dyz.*py+Dzz.*pz;

Opx=Oxy.*py+Oxz.*pz;
Opy=-Oxy.*px+Oyz.*pz;
Opz=-Oxz.*px-Oyz.*py;

sx=epsl*Dpx+Opx;
sy=epsl*Dpy+Opy;
sz=epsl*Dpz+Opz;

fpx=sx+gamap*hx-vx.*dxpx-vy.*dypx-vz.*dzpx;
fpy=sy+gamap*hy-vx.*dxpy-vy.*dypy-vz.*dzpy;
fpz=sz+gamap*hz-vx.*dxpz-vy.*dypz-vz.*dzpz;

end