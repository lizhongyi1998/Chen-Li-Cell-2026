function [rho,vx,vy,vz]=macrovariable(f,Fx,Fy,Fz,const)

rho=zeros(const.Nx,const.Ny,const.Nz);
vx=zeros(const.Nx,const.Ny,const.Nz);
vy=zeros(const.Nx,const.Ny,const.Nz);
vz=zeros(const.Nx,const.Ny,const.Nz);

ex=const.ex;
ey=const.ey;
ez=const.ez;


for i=1:19
    rho=rho+f(:,:,:,i);
    vx=vx+f(:,:,:,i)*ex(i);
    vy=vy+f(:,:,:,i)*ey(i);
    vz=vz+f(:,:,:,i)*ez(i);
end
vx=vx+Fx*const.dt/2;
vy=vy+Fy*const.dt/2;
vz=vz+Fz*const.dt/2;

vx=vx./rho;
vy=vy./rho;
vz=vz./rho;

end