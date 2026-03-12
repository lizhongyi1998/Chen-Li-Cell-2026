clear

const.Nx=30;
const.Ny=30;
const.Nz=50;
const.dx=1;
const.dy=1;
const.dz=1;
const.dt=1;
const.A0 = 0.02;
const.K = 0.02;
const.epsl = 1.2;
const.niu = 2/3;
const.zeta = 0.006;
const.zetac = 0.003;
const.rho = 1;
const.tao = 3*const.niu/const.rho+1/2;
const.gamap = 0.2;
const.Afai=0.01;
const.Kfai=0.04;
const.W=0.04;
const.gamafai=0.05;

R=12;
x0=const.Nx/2;
y0=const.Ny/2;
d1=2;
d2=4;

step = 10000;
Ns = 200;

const.ex=[0 1 -1 0 0 0 0 1 -1 1 -1 0 0 1 -1 1 -1 0 0];
const.ey=[0 0 0 1 -1 0 0 1 -1 0 0 1 -1 -1 1 0 0 1 -1];
const.ez=[0 0 0 0 0 1 -1 0 0 1 -1 1 -1 0 0 -1 1 -1 1];
const.w=[1/3 1/18 1/18 1/18 1/18 1/18 1/18 1/36 1/36 1/36 1/36 1/36 1/36 1/36 1/36 1/36 1/36 1/36 1/36];

[x,y,z] = ndgrid(const.dx*(1:const.Nx),const.dy*(1:const.Ny),const.dz*(1:const.Nz));
fb=(x-x0).^2+(y-y0).^2-R^2;
fb=-fb;
fb(find(fb>=0))=1;
fb(find(fb<0))=0;
reg1=fb;
reg0=1-fb;
reg0(:,:,d1+d2+1:end)=0;
reg0(:,:,1:d1)=1;
reg1(:,:,d1+d2+1:end)=0;
reg1(:,:,1:d1)=0;
const.reg0=reg0;
const.reg1=reg1;
const.W=const.W*(1-reg0-reg1);

theta0=atan2(y-y0,x-x0)+random('normal', pi, 0, const.Nx,const.Ny,const.Nz);
psi0=random('normal', 0.49*pi, 0, const.Nx,const.Ny,const.Nz);
px=sin(psi0).*cos(theta0).*reg1;
py=sin(psi0).*sin(theta0).*reg1;
pz=cos(psi0).*reg1;

fai= ones(const.Nx,const.Ny,const.Nz).*reg1;
vx = zeros(const.Nx,const.Ny,const.Nz);
vy = zeros(const.Nx,const.Ny,const.Nz);
vz = zeros(const.Nx,const.Ny,const.Nz);
rho = const.rho * ones(const.Nx,const.Ny,const.Nz);
f=feq_i(rho,vx,vy,vz,const);

rho_evol = zeros(const.Nx,const.Ny,const.Nz, floor(step / Ns) + 1);
vx_evol = zeros(const.Nx,const.Ny,const.Nz, floor(step / Ns) + 1);
vy_evol = zeros(const.Nx,const.Ny,const.Nz, floor(step / Ns) + 1);
vz_evol = zeros(const.Nx,const.Ny,const.Nz, floor(step / Ns) + 1);
fai_evol = zeros(const.Nx,const.Ny,const.Nz, floor(step / Ns) + 1);
px_evol = zeros(const.Nx,const.Ny,const.Nz, floor(step / Ns) + 1);
py_evol = zeros(const.Nx,const.Ny,const.Nz, floor(step / Ns) + 1);
pz_evol = zeros(const.Nx,const.Ny,const.Nz, floor(step / Ns) + 1);

rho_evol(:, :,:, 1) = rho;
vx_evol(:, :,:, 1) = vx;
vy_evol(:, :,:, 1) = vy;
vz_evol(:, :,:, 1) = vz;
fai_evol(:, :,:, 1) = fai;
px_evol(:, :,:, 1) = px;
py_evol(:, :,:, 1) = py;
pz_evol(:, :,:, 1) = pz;

save_k=1;
for ii=1:step

     [f,rho,vx,vy,vz,fai,px,py,pz]=LBMupdate(f,rho,vx,vy,vz,fai,px,py,pz,const);

     if(mod(ii, Ns) == 0)
        rho_evol(:,:,:, save_k + 1) = rho;
        vx_evol(:,:,:, save_k + 1) = vx;
        vy_evol(:,:,:, save_k + 1) = vy;
        vz_evol(:,:,:, save_k + 1) = vz;
        fai_evol(:,:,:, save_k + 1) = fai;
        px_evol(:,:,:, save_k + 1) = px;
        py_evol(:,:,:, save_k + 1) = py;
        pz_evol(:,:,:, save_k + 1) = pz;
        save_k = save_k + 1;
        save_k
     end

end

save('data.mat','f','rho_evol','vx_evol','vy_evol','vz_evol','fai_evol','px_evol','py_evol','pz_evol','const')