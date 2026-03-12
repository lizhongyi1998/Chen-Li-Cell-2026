function DZ=DiffZ3D(u,dz)

h=zeros(3,3,3);
h(2,2,1)=-0.5/dz;
h(2,2,3)=0.5/dz;
DZ = imfilter(u, h, 'circular');

end