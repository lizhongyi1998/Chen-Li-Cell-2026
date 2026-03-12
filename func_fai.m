function ffai=func_fai(vx,vy,vz,fai,const,varb)

dx=const.dx;
dy=const.dy;
dz=const.dz;
gamafai=const.gamafai;

dxfai=varb.dxfai;
dyfai=varb.dyfai;
dzfai=varb.dzfai;

Dxx=varb.Dxx;
Dyy=varb.Dyy;
Dzz=varb.Dzz;

miu=varb.miu;

ffai=gamafai*Laplacian3D(miu,dx,dy,dz)-vx.*dxfai-vy.*dyfai-vz.*dzfai-(Dxx+Dyy+Dzz).*fai;

end

