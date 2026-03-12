function DX=DiffX3D(u,dx)

h=zeros(3,3,3);
h(1,2,2)=-0.5/dx;
h(3,2,2)=0.5/dx;
DX = imfilter(u, h, 'circular');

end
