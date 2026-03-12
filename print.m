clear

load('data.mat')

ii=40;

fai=fai_evol(:,:,:,ii);
fai1=fai;
fai1(find(fai1>0.5))=1;
fai1(find(fai1<=0.5))=0;

xlim([0,const.Nx])
ylim([0,const.Ny])
zlim([0,const.Nz])
view([1,1,0.6])

xlabel('x')
ylabel('y')
zlabel('z')

% AxesP = [0.1 0.1 0.8 0.8*const.Nz/const.Nx];
AxesP = [0.1 0.1 0.8*const.Nx/const.Nz 0.8];
FigP = [0.25 0.1 0.5 0.8];
set(gca,'position',AxesP);
set(gcf,'unit','normalized','position',FigP,'color','w');

[xM,yM,zM] = ndgrid(const.dx*(1:const.Nx),const.dy*(1:const.Ny),const.dz*(1:const.Nz));

fai(:,:,1:5)=0;
p = patch(isosurface(xM,yM,zM,fai,0.5));
hold on

p.FaceColor = [70,180,250]/255;
p.EdgeColor = 'none';
camlight 
lighting gouraud
alpha(p,0.5)

reg0=const.reg0;
reg0(:,:,1:5)=1;
p = patch(isosurface(xM,yM,zM,reg0,0.5));
hold on

p.FaceColor = [180,180,200]/255;
p.EdgeColor = 'none';
camlight 
lighting gouraud
alpha(p,0.9)

axis off
axis equal


