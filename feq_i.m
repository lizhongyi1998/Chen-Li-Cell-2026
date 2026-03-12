function f=feq_i(rho,vx,vy,vz,const)

f=zeros(const.Nx,const.Ny,const.Nz,19);
ex=const.ex;
ey=const.ey;
ez=const.ez;
w=const.w;

v2=vx.^2+vy.^2+vz.^2;

for i=1:19
    ev=ex(i)*vx+ey(i)*vy+ez(i)*vz;
    f(:,:,:,i)=w(i)*rho.*(1+3*ev+4.5*ev.^2-1.5*v2);
end

end