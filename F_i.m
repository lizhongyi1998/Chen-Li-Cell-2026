function F=F_i(vx,vy,vz,Fx,Fy,Fz,const)

F=zeros(const.Nx,const.Ny,const.Nz,19);

ex=const.ex;
ey=const.ey;
ez=const.ez;
w=const.w;

for i=1:19
    F(:,:,:,i)=(1-1/(2*const.tao))*w(i)*(3*((ex(i)-vx).*Fx+(ey(i)-vy).*Fy+(ez(i)-vz).*Fz)...
        +9*(ex(i)*vx+ey(i)*vy+ez(i)*vz).*(ex(i)*Fx+ey(i)*Fy+ez(i)*Fz));
end

end