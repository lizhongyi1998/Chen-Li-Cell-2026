function DY=DiffY3D(u,dy)

h=zeros(3,3,3);
h(2,1,2)=-0.5/dy;
h(2,3,2)=0.5/dy;
DY = imfilter(u, h, 'circular');

end