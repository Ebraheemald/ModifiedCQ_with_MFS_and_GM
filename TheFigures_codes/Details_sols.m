%=============================================================
% The CQ error as n (number of steps) increasing, 
% where dt = T/N and t = dt*(0:N).
% =============================================================
% T      : Final time
% N      : Number of time steps to compute the soloution 
% Ne    : Number of time steps to compute the exact soloution 
% dt     : Time step for the solution
% dt_e   : Time step for the exact solution
% lag    : Time lag
% X      : Test points
% a2     : Angle of incidence
% a      : a=1 Backward Euler(BE), a=2 the second order backward difference
%formula (BDF2), a=3 Trapezoid rule (TR) for multistep methods.

%====================== sol = 1 ==============================
% The exterior problem for the unit circle
%=============================================================

% M      : Number of collocation points
% Me      : Number of collocation points for the exact soloution
% K      : Number of source points
% Ke     : Number of source points for the exact soloution
% x      : Collocation points
% y      : Source points
% R      : Radius of the circle for the source points

% We have chosen these parameters:
% K = 1000; M = 2000;   Ke = 1500; Me = 3000
% T = 10;  a2 =  -pi/2; lag = 4;   Ne = 4096
% N = [8, 16, 32, 64, 128, 256, 512, 1024]

%====================== sol = 2 ==============================
% Exterior problem for the two Ellipse
%=============================================================

% M      : Number of collocation points
% K      : Number of source points
% Ke     : Number of source points for the exact soloution
% x      : Collocation points
% y      : Source points
% R      : Radius of the circle for the source points

% We have chosen these parameters:
% K = 2000; M = 4000;    Ke = 3000; Me = 6000
% T = 10;  a2 =  -pi/2; lag = 4;    Ne = 4096
% N = [8, 16, 32, 64, 128, 256, 512, 1024]

%====================== sol = 3 ==============================
% The exterior problem for the semicircles
%=============================================================

% M      : boundary points
% Me     : boundary points for the exact soloution

% We have chosen these parameters:
%  M = 1000;   Me = 1500; T = 10;
% a2 =  pi/4; lag = 4;   Ne = 4096
%  N = [8, 16, 32, 64, 128, 256, 512, 1024]

%====================== sol = 4 ==============================
% The exterior problem for the semicircles
%=============================================================

% M      : boundary points
% Me     : boundary points for the exact soloution

% We have chosen these parameters:
%  M = 2000;   Me = 3000; T = 10;
% a2 =  pi/4; lag = 4;   Ne = 4096
%  N = [8, 16, 32, 64, 128, 256, 512, 1024]
