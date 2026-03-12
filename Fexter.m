function [Fx,Fy,Fz]=Fexter(fai,px,py,pz,const,varb)

dx=const.dx;
dy=const.dy;
dz=const.dz;
K=const.K;
epsl=const.epsl;
zeta=const.zeta;
zetac=const.zetac;
Kfai=const.Kfai;

dxfai=varb.dxfai;
dyfai=varb.dyfai;
dzfai=varb.dzfai;

dpdpxx=varb.dpdpxx;
dpdpyy=varb.dpdpyy;
dpdpzz=varb.dpdpzz;
dpdpxy=varb.dpdpxy;
dpdpxz=varb.dpdpxz;
dpdpyz=varb.dpdpyz;

hx=varb.hx;
hy=varb.hy;
hz=varb.hz;

sgmsxx=-epsl*hx.*px;
sgmsyy=-epsl*hy.*py;
sgmszz=-epsl*hz.*pz;
sgmsxy=-epsl/2*(hx.*py+hy.*px);
sgmsyz=-epsl/2*(hy.*pz+hz.*py);
sgmsxz=-epsl/2*(hx.*pz+hz.*px);

Trsgms=sgmsxx+sgmsyy+sgmszz;
sgmsxx=sgmsxx-Trsgms/3;
sgmsyy=sgmsyy-Trsgms/3;
sgmszz=sgmszz-Trsgms/3;

sgmaxy=1/2*(px.*hy-py.*hx);
sgmayz=1/2*(py.*hz-pz.*hy);
sgmaxz=1/2*(px.*hz-pz.*hx);

sgmexx=-K*dpdpxx;
sgmeyy=-K*dpdpyy;
sgmezz=-K*dpdpzz;
sgmexy=-K*dpdpxy;
sgmeyz=-K*dpdpyz;
sgmexz=-K*dpdpxz;

sgmactxx=-zeta*px.*px;
sgmactyy=-zeta*py.*py;
sgmactzz=-zeta*pz.*pz;
sgmactxy=-zeta*px.*py;
sgmactyz=-zeta*py.*pz;
sgmactxz=-zeta*px.*pz;

Trsgmact=sgmactxx+sgmactyy+sgmactzz;
sgmactxx=sgmactxx-Trsgmact/3;
sgmactyy=sgmactyy-Trsgmact/3;
sgmactzz=sgmactzz-Trsgmact/3;

Cppxx=varb.Cppxx;
Cppyy=varb.Cppyy;
Cppzz=varb.Cppzz;
Cppxy=varb.Cppxy;
Cppyx=varb.Cppyx;
Cppyz=varb.Cppyz;
Cppzy=varb.Cppzy;
Cppxz=varb.Cppxz;
Cppzx=varb.Cppzx;

sgmacxx=-zetac*Cppxx;
sgmacyy=-zetac*Cppyy;
sgmaczz=-zetac*Cppzz;
sgmacxy=-zetac*Cppxy;
sgmacyx=-zetac*Cppyx;
sgmacyz=-zetac*Cppyz;
sgmaczy=-zetac*Cppzy;
sgmacxz=-zetac*Cppxz;
sgmaczx=-zetac*Cppzx;

sgmcxx=-Kfai*dxfai.*dxfai;
sgmcyy=-Kfai*dyfai.*dyfai;
sgmczz=-Kfai*dzfai.*dzfai;
sgmcxy=-Kfai*dxfai.*dyfai;
sgmcxz=-Kfai*dxfai.*dzfai;
sgmcyx=-Kfai*dxfai.*dyfai;
sgmcyz=-Kfai*dzfai.*dyfai;
sgmczx=-Kfai*dxfai.*dzfai;
sgmczy=-Kfai*dzfai.*dyfai;

Trsgmc=sgmcxx+sgmcyy+sgmczz;
sgmcxx=sgmcxx-Trsgmc/3;
sgmcyy=sgmcyy-Trsgmc/3;
sgmczz=sgmczz-Trsgmc/3;

sgmxx=sgmsxx+sgmexx+sgmactxx+sgmacxx;
sgmyy=sgmsyy+sgmeyy+sgmactyy+sgmacyy;
sgmzz=sgmszz+sgmezz+sgmactzz+sgmaczz;
sgmxy=sgmsxy+sgmaxy+sgmexy+sgmactxy+sgmacxy;
sgmyx=sgmsxy-sgmaxy+sgmexy+sgmactxy+sgmacyx;
sgmyz=sgmsyz+sgmayz+sgmeyz+sgmactyz+sgmacyz;
sgmzy=sgmsyz-sgmayz+sgmeyz+sgmactyz+sgmaczy;
sgmxz=sgmsxz+sgmaxz+sgmexz+sgmactxz+sgmacxz;
sgmzx=sgmsxz-sgmaxz+sgmexz+sgmactxz+sgmaczx;

dxsgmxx=DiffX3D(sgmxx,dx);
dxsgmyx=DiffX3D(sgmyx,dx);
dxsgmzx=DiffX3D(sgmzx,dx);
dysgmxy=DiffY3D(sgmxy,dy);
dysgmyy=DiffY3D(sgmyy,dy);
dysgmzy=DiffY3D(sgmzy,dy);
dzsgmxz=DiffZ3D(sgmxz,dz);
dzsgmyz=DiffZ3D(sgmyz,dz);
dzsgmzz=DiffZ3D(sgmzz,dz);

Fx=dxsgmxx+dysgmxy+dzsgmxz;
Fy=dxsgmyx+dysgmyy+dzsgmyz;
Fz=dxsgmzx+dysgmzy+dzsgmzz;

dxsgmcxx=DiffX3D(sgmcxx,dx);
dxsgmcyx=DiffX3D(sgmcyx,dx);
dxsgmczx=DiffX3D(sgmczx,dx);
dysgmcxy=DiffY3D(sgmcxy,dy);
dysgmcyy=DiffY3D(sgmcyy,dy);
dysgmczy=DiffY3D(sgmczy,dy);
dzsgmcxz=DiffZ3D(sgmcxz,dz);
dzsgmcyz=DiffZ3D(sgmcyz,dz);
dzsgmczz=DiffZ3D(sgmczz,dz);

Fcx=dxsgmcxx+dysgmcxy+dzsgmcxz;
Fcy=dxsgmcyx+dysgmcyy+dzsgmcyz;
Fcz=dxsgmczx+dysgmczy+dzsgmczz;

Fcx=Fcx.*(1-const.reg0-const.reg1);
Fcy=Fcy.*(1-const.reg0-const.reg1);
Fcz=Fcz.*(1-const.reg0-const.reg1);

Fx=Fx+Fcx;
Fy=Fy+Fcy;
Fz=Fz+Fcz;

end