function varb=func_varb(fai,px,py,pz,const)

dx=const.dx;
dy=const.dy;
dz=const.dz;
A0=const.A0;
K=const.K;
W=const.W;
Afai=const.Afai;
Kfai=const.Kfai;

dxfai=DiffX3D(fai,dx);
dyfai=DiffY3D(fai,dy);
dzfai=DiffZ3D(fai,dz);

dxpx=DiffX3D(px,dx);
dxpy=DiffX3D(py,dx);
dxpz=DiffX3D(pz,dx);
dypx=DiffY3D(px,dy);
dypy=DiffY3D(py,dy);
dypz=DiffY3D(pz,dy);
dzpx=DiffZ3D(px,dz);
dzpy=DiffZ3D(py,dz);
dzpz=DiffZ3D(pz,dz);

dyppxx=2*dypx.*px;
dzppxx=2*dzpx.*px;
dxppyy=2*dxpy.*py;
dzppyy=2*dzpy.*py;
dxppzz=2*dxpz.*pz;
dyppzz=2*dypz.*pz;
dxppxy=dxpx.*py+px.*dxpy;
dyppxy=dypx.*py+px.*dypy;
dzppxy=dzpx.*py+px.*dzpy;
dxppxz=dxpx.*pz+px.*dxpz;
dyppxz=dypx.*pz+px.*dypz;
dzppxz=dzpx.*pz+px.*dzpz;
dxppyz=dxpy.*pz+py.*dxpz;
dyppyz=dypy.*pz+py.*dypz;
dzppyz=dzpy.*pz+py.*dzpz;

Cppxx=dyppxz-dzppxy;
Cppyy=dzppxy-dxppyz;
Cppzz=dxppyz-dyppxz;
Cppxy=dyppyz-dzppyy;
Cppyx=dzppxx-dxppxz;
Cppxz=dyppzz-dzppyz;
Cppzx=dxppxy-dyppxx;
Cppyz=dzppxz-dxppzz;
Cppzy=dxppyy-dyppxy;

dpdpxx=dxpx.^2+dxpy.^2+dxpz.^2;
dpdpyy=dypx.^2+dypy.^2+dypz.^2;
dpdpzz=dzpx.^2+dzpy.^2+dzpz.^2;
dpdpxy=dxpx.*dypx+dxpy.*dypy+dxpz.*dypz;
dpdpyz=dypx.*dzpx+dypy.*dzpy+dypz.*dzpz;
dpdpxz=dxpx.*dzpx+dxpy.*dzpy+dxpz.*dzpz;

pdf=px.*dxfai+py.*dyfai+pz.*dzfai;

p2=px.^2+py.^2+pz.^2;
hx=-A0*(-fai+p2).*px+K*Laplacian3D(px,dx,dy,dz);
hy=-A0*(-fai+p2).*py+K*Laplacian3D(py,dx,dy,dz);
hz=-A0*(-fai+p2).*pz+K*Laplacian3D(pz,dx,dy,dz);

hsux=-2*W.*pdf.*dxfai;
hsuy=-2*W.*pdf.*dyfai;
hsuz=-2*W.*pdf.*dzfai;

Wb=-0.2;
Wb=Wb*(const.reg0+const.reg1);
hsubx=-Wb.*dxfai;
hsuby=-Wb.*dyfai;
hsubz=0;

hx=hx+hsux+hsubx;
hy=hy+hsuy+hsuby;
hz=hz+hsuz+hsubz;

miu=Afai*(fai-fai.^2).*(1-2*fai)-A0/2*p2-Kfai*Laplacian3D(fai,dx,dy,dz);

varb.dxfai=dxfai;
varb.dyfai=dyfai;
varb.dzfai=dzfai;
varb.dxpx=dxpx;
varb.dxpy=dxpy;
varb.dxpz=dxpz;
varb.dypx=dypx;
varb.dypy=dypy;
varb.dypz=dypz;
varb.dzpx=dzpx;
varb.dzpy=dzpy;
varb.dzpz=dzpz;
varb.Cppxx=Cppxx;
varb.Cppyy=Cppyy;
varb.Cppzz=Cppzz;
varb.Cppxy=Cppxy;
varb.Cppyx=Cppyx;
varb.Cppyz=Cppyz;
varb.Cppzy=Cppzy;
varb.Cppxz=Cppxz;
varb.Cppzx=Cppzx;
varb.dpdpxx=dpdpxx;
varb.dpdpyy=dpdpyy;
varb.dpdpzz=dpdpzz;
varb.dpdpxy=dpdpxy;
varb.dpdpxz=dpdpxz;
varb.dpdpyz=dpdpyz;
varb.pdf=pdf;
varb.hx=hx;
varb.hy=hy;
varb.hz=hz;
varb.miu=miu;

end