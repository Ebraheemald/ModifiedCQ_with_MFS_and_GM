%%
% The exterior problem for the nonconvex domain (semicircles)

M = 2000;

%== connected semicircles ==============
p=semicircles_f(M);

% bounding box containing the semcircles
b = [-4-4i 4-4i 4+4i -4+4i];
pgon = polyshape({real(b), real(p)},{imag(b),imag(p)}); %singly connected
plot(pgon)

%%
tr = triangulation(pgon);
model = createpde;
tnodes = tr.Points';
telements = tr.ConnectivityList';

geometryFromMesh(model,tnodes,telements);
generateMesh(model,'Hmax',0.05);
pdemesh(model);

X = model.Mesh.Nodes(1,:).';
Y = model.Mesh.Nodes(2,:).';
tri = model.Mesh.Elements(1:3,:).';
%

load('GM_Semicircles_data_w5_pi4')
N = N1(end);
ts = (0:N)*T/N;
load('GM_Semicircles_alpha2_w5_BDF2_pi4_z2.mat')
% z1=CQ_eval_sys_mod(N,T,(X+1i*Y).',alpha,2,p);
% z2=z1(:,t2), where t1=2^8:2^7:N-2^7, and t2=t1+1
z2=real(z2);

close all
time=zeros(length(t1),1); count=0;
for j = t1
    figure; count=count+1;
    time(count)=(j/N)*T;
    Z = real(z2(:,count)-g((X+1i*Y),ts(j+1)));
    trisurf(tri,X,Y,Z); axis equal;
    caxis([-.25,.25])
    view(2)
    grid off
    shading interp
    pause(.1)
    xlim([-4 4])
    ylim([-4 4])
end