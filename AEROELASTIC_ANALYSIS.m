clear
close all   
clc
%%
% Data
lamda = 34; % Sweep angle of the wing
b = 11; % Semi span (m)
bbar = b/cosd(lamda);
c1 = 4.5; % Root cord length (m)
c2 = 2.5; % Tip cord length (m)
Mach = 0.7; % Mach number
Cla0 = 2*pi; % Aerodynamic coefficients
Cla = Cla0/(sqrt(1-Mach^2));
Clabar = Cla/cosd(lamda);
E = 0.25;
zeta = 0.02; % Structural Damping

% Symbolic Matlab
syms ybar a rho Vg

% Data EI, GJ, m, Inertia wrt EA (It), xcg
EI = 10^6*(6.7 - 5*(ybar/bbar));
GJ = 10^6*(9.38 - 6.25*(ybar/bbar))*0.5;
m = 68 - 23*(ybar/bbar);
It = 160 - 51*(ybar/bbar);
xcg = 0.38 - 0.5*(ybar/bbar);

% Shape Function
Nw = [(ybar/bbar)^2 (ybar/bbar)^3 (ybar/bbar)^4 (ybar/bbar)^5 (ybar/bbar)^6];
Nt = [sin(pi/2*(ybar/bbar)) sin(pi*(ybar/bbar)) sin(3/2*pi*(ybar/bbar)) sin(2*pi*(ybar/bbar))];

% First Derivative of Shape Function
dNw = diff(Nw,ybar);
dNt = diff(Nt,ybar);

% Second Derivative of Shape Function
ddNw = diff(dNw,ybar);
ddNt = diff(dNt,ybar);

% Stiffness Matrix (STRUCTURE)
Kww = int(ddNw'*EI*ddNw,ybar,[0 bbar]);
Ktt = int(dNt'*GJ*dNt,ybar,[0 bbar]);

% Mass Matrix
Mww = int(Nw'*m*Nw,ybar,[0 bbar]);
Mwt = int(Nw'*m*xcg*Nt,ybar,[0 bbar]);
Mtw = int(Nt'*m*xcg*Nw,ybar,[0 bbar]);
Mtt = int(Nt'*It*Nt,ybar,[0 bbar]);

% Quasi-Steady Aerodynamics Approach
Uinf = Mach*a;
q = 0.5*rho*Uinf^2;
qbar = q*cosd(lamda)^2;
c = c1 - (c1-c2)*(ybar/bbar);
cbar = c*cosd(lamda);
e = (0.43-0.25)*c;
ebar = e*cosd(lamda);
Kwwa = -int(Nw'*qbar*cbar*Clabar*dNw*tand(lamda),ybar,[0 bbar]) + int(dNw'*tand(lamda)*ebar*qbar*cbar*Clabar*dNw*tand(lamda),ybar,[0 bbar]);
Kwta = int(Nw'*qbar*cbar*Clabar*Nt,ybar,[0 bbar]) - int(dNw'*tand(lamda)*ebar*qbar*cbar*Clabar*Nt,ybar,[0 bbar]);
Ktwa = -int(Nt'*ebar*qbar*cbar*Clabar*dNw*tand(lamda),ybar,[0 bbar]);
Ktta = int(Nt'*ebar*qbar*cbar*Clabar*Nt,ybar,[0 bbar]);
Cwwa = -int(Nw'*qbar*cbar*Clabar*Nw*(1/Uinf)*(1/cosd(lamda)),ybar,[0 bbar]) + int(Nw'*qbar*cbar*Clabar*(ebar-cbar/2)*dNw*(1/Uinf)*(sind(lamda)/cosd(lamda)^2),ybar,[0 bbar]) + int(dNw'*ebar*qbar*cbar*Clabar*Nw*(1/Uinf)*(sind(lamda)/cosd(lamda)^2),ybar,[0 bbar]) - int(dNw'*ebar*qbar*cbar*Clabar*(ebar-cbar/2)*dNw*(1/Uinf)*(sind(lamda)^2/cosd(lamda)^3),ybar,[0 bbar]);
Cwta = -int(Nw'*qbar*cbar*Clabar*(ebar-cbar/2)*Nt*(1/Uinf)*(1/cosd(lamda)),ybar,[0 bbar]) + int(dNw'*ebar*qbar*cbar*Clabar*(ebar-cbar/2)*Nt*(1/Uinf)*(sind(lamda)/cosd(lamda)^2),ybar,[0 bbar]);
Ctwa = -int(Nt'*ebar*qbar*cbar*Clabar*Nw*(1/Uinf)*(1/cosd(lamda)),ybar,[0 bbar]) + int(Nt'*ebar*qbar*cbar*Clabar*(ebar-cbar/2)*dNw*(1/Uinf)*(sind(lamda)/cosd(lamda)^2),ybar,[0 bbar]);
Ctta = -int(Nt'*ebar*qbar*cbar*Clabar*(ebar-cbar/2)*Nt*(1/Uinf)*(1/cosd(lamda)),ybar,[0 bbar]);

% Forcing Term by Gust
Fw = -int(Nw'*qbar*cbar*Clabar/cosd(lamda),ybar,[0 bbar]) + int(dNw'*ebar*qbar*cbar*Clabar*tand(lamda)/cosd(lamda),ybar,[0 bbar]);
Ft = -int(Nt'*ebar*qbar*cbar*Clabar/cosd(lamda),ybar,[0 bbar]);

%% Defining the matrices
[~,nw] = size(Nw);
[~,nt] = size(Nt);
M = double([Mww -Mwt;-Mtw Mtt]);
K = double([Kww zeros(nw,nt);zeros(nt,nw) Ktt]);
Ka = vpa([Kwwa Kwta;Ktwa Ktta]);
Ca = vpa([Cwwa Cwta;Ctwa Ctta]);
F = vpa([Fw;Ft]);

% Check the frequency less than 10 Hz
[V,D] = eig(K,M);
eigenval = sqrt(diag(D));
eigenval = eigenval/(2*pi); % convert eigenvalues in Hz
freq = 10;
nfreq = 0;
while eigenval(nfreq+1) < 10
    nfreq = nfreq+1;
end

fprintf('Number of Eigenvalues less than 10 Hz is %.f\n\n', nfreq)

% System matrices using the frequencies below 10Hz
MM = V(:,1:nfreq)'*M*V(:,1:nfreq);
KK = V(:,1:nfreq)'*K*V(:,1:nfreq);
KKa = V(:,1:nfreq)'*Ka*V(:,1:nfreq);
CC = 2*zeta*MM*sqrt(D(1:nfreq,1:nfreq));
CCa = V(:,1:nfreq)'*Ca*V(:,1:nfreq);
FF = V(:,1:nfreq)'*F;


%% Finding the Velocity and Dynamics Pressure at Flutter Condition

A = vpa([-inv(MM)*(CC - CCa) -inv(MM)*(KK - KKa);eye(nfreq) zeros(nfreq)]);
h = linspace(0,40000,2000);
[T, aisa, P, rhoisa] = atmoscoesa(h);
DDf = zeros(length(A),length(h));
qqftest = zeros(1,length(h));

for i = 1:length(h)
    AA = double(subs(A,{a rho},{aisa(i) rhoisa(i)}));
    qqftest(i) = 0.5*rhoisa(i)*(aisa(i)*Mach)^2;
    [Vf, Df] = eig(AA);
    DDf(:,i) = diag(Df);
end
i = 1;
while max(real(DDf(:,i))) > 0
    i = i+1;
end

id = i;
Uf = aisa(id)*Mach;
qf = 0.5*rhoisa(id)*Uf^2;

fprintf('Velocity at flutter condition in %.2f m/s\n\n', Uf)
fprintf('Dynamics pressure at flutter condition in %.2f Pa\n\n', qf)

%% Bending at root due to 1-cosine gust
d = 16;
qqf = (1 - d/100)*qf;
err = zeros(1,length(i));

for i = 1:length(h)
    err(i) = abs(qqftest(i) - qqf);
end

[error,idd] = min(err);

UUf = aisa(idd)*Mach;
CCa = subs(CCa,{a rho}, {aisa(idd) rhoisa(idd)});
KKa = subs(KKa,{a rho}, {aisa(idd) rhoisa(idd)});
FF = subs(FF,{a rho}, {aisa(idd) rhoisa(idd)});


odefun = @(t,x) Gust(t, x, nfreq, MM, CC, CCa, KK, KKa, FF, UUf);

wg0 = 6.35;
Lg = 15.24;

tspan = 0:0.001:5;
tspan_gust = 0:0.001:Lg/UUf;
t_gust = [tspan_gust zeros(1,length(tspan)-length(tspan_gust))];

[t,z] = ode45(odefun,tspan,zeros(2*nfreq,1));

ddz = double([-inv(MM)*(CC-CCa) -inv(MM)*(KK-KKa)])*z' + double(inv(MM)\FF)*(wg0/2)*(1-cos(2*pi*UUf*t_gust/Lg))*(1/UUf);

qw = V(1:nw,1:nfreq)*z(:,nfreq+1:end)';

bending = double(subs(EI,ybar,0)*subs(ddNw,ybar,0)*qw);
bending2 = double(subs(EI,ybar,0)*subs(ddNw,ybar,0)*qw);


plot(t,bending); grid on;
legend('Acceleration Mode')
fprintf('Maximum bending at root using direct recovery %.2f Nm\n\n', max(abs(bending)))

%%
function  dx = Gust(t, x, nfreq, MM, CC, CCa, KK, KKa, FF, UUf)

   wg0 = 6.35;
   Lg = 15.24;

   Ass = double([-inv(MM)*(CC-CCa) -inv(MM)*(KK-KKa);eye(nfreq) zeros(nfreq)]);
   Bss = double([inv(MM)*FF; zeros(nfreq,1)]);
   if t < Lg/UUf
       u = (wg0/2)*(1-cos(2*pi*UUf*t/Lg))/UUf;
   else 
       u = 0;
   end
   dx = Ass*x + Bss*u;

end