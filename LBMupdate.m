function [f,rho,vx,vy,vz,fai,px,py,pz]=LBMupdate(f0,rho0,vx0,vy0,vz0,fai0,px0,py0,pz0,const)

dt=const.dt;
tao=const.tao;
reg0=const.reg0;
reg1=const.reg1;

varb=func_varb(fai0,px0,py0,pz0,const);
[Fx0,Fy0,Fz0]=Fexter(fai0,px0,py0,pz0,const,varb);
[rho_t0,vx_t0,vy_t0,vz_t0]=macrovariable(f0,Fx0,Fy0,Fz0,const);
[fpx1,fpy1,fpz1,varb]=func_p(vx_t0,vy_t0,vz_t0,px0,py0,pz0,const,varb);
ffai1=func_fai(vx_t0,vy_t0,vz_t0,fai0,const,varb);
Fi0=F_i(vx_t0,vy_t0,vz_t0,Fx0,Fy0,Fz0,const);

pNx1_try = px0 + 0.5 * fpx1*dt;
pNy1_try = py0 + 0.5 * fpy1*dt;
pNz1_try = pz0 + 0.5 * fpz1*dt;
px1_try = pNx1_try + 0.5 * fpx1*dt;
py1_try = pNy1_try + 0.5 * fpy1*dt;
pz1_try = pNz1_try + 0.5 * fpz1*dt;

faiN1_try=fai0+0.5*ffai1*dt;
fai1_try=faiN1_try+0.5*ffai1*dt;

faiN1_try=faiN1_try.*(1-reg0);
faiN1_try=faiN1_try.*(1-reg1)+reg1;
fai1_try=fai1_try.*(1-reg0);
fai1_try=fai1_try.*(1-reg1)+reg1;

fi_eq = feq_i(rho_t0,vx_t0,vy_t0,vz_t0,const);
fi_n = f0 + 0.5 * (1/tao * (fi_eq - f0) + Fi0)*dt;
fi_f = fi_n + 0.5 * (1/tao * (fi_eq - f0) + Fi0)*dt;

fi_fstreaming = stream(fi_f,const);
fi_nstreaming = stream(fi_n,const);

varb=func_varb(fai1_try,px1_try,py1_try,pz1_try,const);
[Fx1,Fy1,Fz1]=Fexter(fai1_try,px1_try,py1_try,pz1_try,const,varb);
[rho_t1_try, vx_t1_try, vy_t1_try, vz_t1_try] =macrovariable(fi_fstreaming,Fx1,Fy1,Fz1,const);
[fpx2,fpy2,fpz2,varb]=func_p(vx_t1_try,vy_t1_try,vz_t1_try,px1_try,py1_try,pz1_try,const,varb);
ffai2=func_fai(vx_t1_try,vy_t1_try,vz_t1_try,fai1_try,const,varb);
Fi1=F_i(vx_t1_try, vy_t1_try, vz_t1_try,Fx1,Fy1,Fz1,const);

px1 = pNx1_try + 0.5 * fpx2*dt;
py1 = pNy1_try + 0.5 * fpy2*dt;
pz1 = pNz1_try + 0.5 * fpz2*dt;

fai1=faiN1_try+0.5*ffai2*dt;

fai1=fai1.*(1-reg0);
fai1=fai1.*(1-reg1)+reg1;

fi_eq = feq_i(rho_t1_try, vx_t1_try, vy_t1_try, vz_t1_try,const);
fi_f = fi_nstreaming + 0.5 * (1/tao* (fi_eq - fi_fstreaming) + Fi1)*dt;

rho=rho_t1_try;
vx=vx_t1_try;
vy=vy_t1_try;
vz=vz_t1_try;
f=fi_f;
fai=fai1;
px=px1;
py=py1;
pz=pz1;

end