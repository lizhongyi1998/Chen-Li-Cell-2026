function del2u = Laplacian3D(u, dx, dy,dz)

h=zeros(3,3,3);
h(2,2,2)=-2*(1/dx^2+1/dy^2+1/dz^2);
h(1,2,2)=1/dx^2;
h(3,2,2)=1/dx^2;
h(2,1,2)=1/dy^2;
h(2,3,2)=1/dy^2;
h(2,2,1)=1/dz^2;
h(2,2,3)=1/dz^2;
del2u = imfilter(u, h, 'circular');

end