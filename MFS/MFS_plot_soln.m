%%
% The exterior problem for the two ellipses

M = 2000;
z = exp(1i*2*pi*(0:M-1)/M);
x = (z+.2./z)*exp(1i*pi/2); x = x*.5; x=x.';
x1 = x-2; x2 = x+2;

% bounding box containing the two ellipses
b = [-4-4i 4-4i 4+4i -4+4i];
pgon = polyshape({real(b), real(x1), real(x2)},{imag(b),imag(x1), imag(x2)}); %doubly connected
plot(pgon)

%%
tr = triangulation(pgon);
model = createpde;
tnodes = tr.Points';
telements = tr.ConnectivityList';

geometryFromMesh(model,tnodes,telements);
generateMesh(model,'Hmax',.05);
pdemesh(model);

X = model.Mesh.Nodes(1,:).';
Y = model.Mesh.Nodes(2,:).';
tri = model.Mesh.Elements(1:3,:).';
%

load('TwoEllipses_BDF2_alpha2_w5')
N = N1(end);
ts = (0:N)*T/N;
K = K/2;
ys1 = R*exp(1i*2*pi*(0:K-1)/K);
ys = (ys1+.2./ys1)*exp(1i*pi/2); y = ys*.5;
y1 = y-2; y2 = y+2; y=[y1 y2]; K=2*K;
load('TwoEllipses_BDF2_alpha2_w5_z2.mat')
% z1=eval_CQ_sys_mod(N,T,X+1i*Y,y,alpha,2);
% z2=z1(:,t2), where t1=2^8:2^7:N-2^7, and t2=t1+1

close all
time=zeros(length(t1),1); count=0;
for j = t1      
    figure; count=count+1;
    time(count)=(j/N)*T;
    Z = z2(:,count)-g(X+1i*Y,ts(j+1)); 
    trisurf(tri,X,Y,Z); axis equal;
    caxis([-.25,.25])
    view(2)
    grid off
    shading interp
    pause(.2)
    xlim([-4 4])
    ylim([-4 4])
end