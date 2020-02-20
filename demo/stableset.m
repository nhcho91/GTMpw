% Demonstration of a region-of-attraction estimation
% for the polynomial longitudinal dynamics.
%
%% Dependencies
%
% Requires
% * pwroaest -- https://github.com/pwpfit/pwroaest
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2018-11-28
% * Changed:    2018-11-28
%
%% Variables, constants, and their units
%
% * |alpha|    :  angle of attack,                                  rad
% * |gamma|    :  flight-path angle,                                rad
% * |gammadot| :  change in flight-path angle,                      rad/s
% * |eta|      :  elevator deflection,                              rad
% * |M|        :  pitch moment, body y-axis,                        lbf-ft
% * |T|        :  thrust,                                           lbf
% * |V|        :  airspeed,                                         ft/s
% * |Vdot|     :  change in airspeed,                               ft/s^2
%
%%

warning('off', 'MATLAB:nargchk:deprecated')

%% Equations of motion
import aero.GtmPiecewise
import eom.GtmLong
import eom.GtmLpoly

% piecewise aerodynamic coefficients 
pw = GtmPiecewise;

% pre-stall boundary
alphapw = pw.alpha0;

% non-linear equations of motion
eom.nl = GtmLong(pw);

% equations of motion
f = @(EOM) @(x,u,varargin) EOM.f(x,u);

% level-flight trim condition (gamma = 0)
import aerootools.findtrim
x0 = [35; 0; 0; 0.1];
[x0,u0] = findtrim(f(eom.nl), x0, zeros(2,1), @(x,u) x(1:2)-x0(1:2));

% polynomial equations of motion
eom.p1 = GtmLpoly(x0,pw.pre);
eom.p2 = GtmLpoly(x0,pw.post);

f1 = f(eom.p1);
f2 = f(eom.p2);

% scaling
% x_sc = D*x
Dx = [
        20              % m/s
        deg2rad(20)     % rad
        deg2rad(150)    % rad
        deg2rad(20)     % rad
].^-1;

% pitch rate damping
% from CSB+2011
Kq = 0.0698;

% elevator constraints
ebnd = deg2rad([-30  20]);

% constraint function
c = @(x,u) (u(1) - sum(ebnd)/2).^2 - (max(ebnd) - sum(ebnd)/2)^2;


%% Polynomials
pvar V alpha gamma q
X = double(eom.nl.X(V,gamma,q,alpha));

% control law
% if LQR feedback is employed, replace by
% K = -lqr(...);    % lqr returns inverted feedback
% U = -K*X;     
U = [Kq;0]*q; 

% polynomial vector fields
F = f1(X+x0, U+u0);     % low-angle of attack
G = f2(X+x0, U+u0);     % high-angle of attack

% polynomial constraint
C = c(X+x0,U+u0);

% clean polynomial
% remove coefficients < 10^-6
% remove terms of order > 6
Fb = cleanpoly(F, 1e-6, 0:5);
Gb = cleanpoly(G, 1e-6, 0:5);

% scaling 
D = diag(Dx);
% Y = D*X
Y  = X;
% ydot = D*xdot = D*f(X) = D*f(Y/D)
Fl = subs(D*Fb, X, D^-1*Y);
Gl = subs(D*Gb, X, D^-1*Y);

% scaled boundary condition
Xi  = alpha;
Phi = subs(alpha + x0(4) - alphapw, X, D^-1*Y);

% input scaling
d = deg2rad(20)^-1;
% scaled constraint
Cl = subs(d^2*C, X, D^-1*Y);

% scaled control
% uncomment if needed
Kl = subs(d*U, X, D^-1*Y);


%% Piecewise safe-set estimation
zV = monomials(Y,2:4);
z1 = monomials(Y,0:1);
z2 = monomials(Y,1:2);
zg = monomials(Y,0:1);
l1 = 1e-6*(Y'*Y);

zi = monomials(Y,0:2);

% options for piecewise ROA estimation
% see documentation for full list of options
ropt = pwroaoptions(Fl, Gl, Phi, Y, 'c', Cl, 'zV', zV, 'p', Y'*Y*1e2, 'betamax', 1000, 'zg', zg, 'zi', zi, 'z2', z2, 'z1', z1, 'L1', l1, 'display', 'on');
% set number of iterations
ropt.NstepBis = 50;

[b,Lya,g,~,it] = pwroaest(ropt);


%% Plot
La = subs(Lya{1},   Y, D*X); La = subs(La,  X(2:end), deg2rad(1)*X(2:end));
Lb = subs(Lya{end}, Y, D*X); Lb = subs(Lb,  X(2:end), deg2rad(1)*X(2:end));
P  = subs(ropt.p,   Y, D*X); P  = subs(P,   X(2:end), deg2rad(1)*X(2:end));
H  = subs(Phi+1e-12*sum(Y), Y, D*X); 
                             H  = subs(H,   X(2:end), deg2rad(1)*X(2:end));
k  = subs(rad2deg(1)*(d^-1*Kl(1)+u0(1))+1e-12*sum(Y), Y, D*X); 
                             k  = subs(k,   X(2:end), deg2rad(1)*X(2:end));

% phugoid plane
L2a = subs(La, X(3:4), zeros(2,1));
L2b = subs(Lb, X(3:4), zeros(2,1));
P2  = subs(P,  X(3:4), zeros(2,1));
k2  = subs(k,  X(3:4), zeros(2,1));

% short-period plane
L1a = subs(La, X(1:2), zeros(2,1));
L1b = subs(Lb, X(1:2), zeros(2,1));
P1  = subs(P,  X(1:2), zeros(2,1));
H1  = subs(H,  X(1:2), zeros(2,1));
k1  = subs(k,  X(1:2), zeros(2,1));

L = round(rad2deg(linspace(ebnd(1),ebnd(2),11)));

figure(1)
clf

subplot(1,2,1)
title('phugoid motion')
hold on
pcontour(L2a, g, [-17 17 -35 35], 'r');
pcontour(L2b, g, [-17 17 -35 35], 'r--');
pcontour(P2,  b, [-17 17 -35 35], 'b');
grid on
[C,h] = ...
pcontour(k2,  L, [-17 17 -35 35], '-.');
clabel(C,h);


subplot(1,2,2)
title('short-period motion')
hold on
pcontour(L1a, g, [-30    9.4 -225 225], 'r');
pcontour(L1b, g, [  9.4 30   -225 225], 'r--');
pcontour(P1,  b, [-30   30   -225 225], 'b');
grid on
pcontour(H1,  0, [-30   30   -225 225], 'k:');

[C,h] = ...
pcontour(k1,  L, [-30   30   -225 225], '-.');
clabel(C,h);

legend('invariant set', 'invariant set (high alpha)', 'ellispoidal shape', 'boundary condition', 'control feedback (deg)');
