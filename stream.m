function f=stream(f0,const)

ex=const.ex;
ey=const.ey;
ez=const.ez;

for i=1:19
    f(:,:,:,i)=circshift(f0(:,:,:,i),[ex(i),ey(i),ez(i),0]);
end

end