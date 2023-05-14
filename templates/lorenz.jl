using PyPlot;
using ODE;

function f(t,r)
    (x,y,z)=r
    dx_dt=σ*(y-x)
    dy_dt=x*(ρ-z)-y
    dz_dt=x*y-β*z
    [dx_dt;dy_dt;dz_dt]
end;

dt=0.001
tf=100.0
t=collect(0:dt:tf)
r0=[0.1;0.0;0.0]
σ=10.0
ρ=28.0
β=8.0/3.0;
(t,pos)=ode23(f,r0,t)
x=map(v -> v[1], pos)
y=map(v -> v[2], pos)
z=map(v -> v[3], pos);
plot3D(x,y,z);
show();
