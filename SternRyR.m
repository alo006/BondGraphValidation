function SternRyR

% Variables

% y(1) = RyRr0;
% y(2) = RyRo0;
% y(3) = RyRi0;
% y(4) = Cai;
% y(5) = Casr;

% Clamp parameters

cajclamp = 8.597401e-5;
casrclamp = 0.1e-1;

% RyR parameters (Shannon/Grandi 2011)

koCa = 10+20*AF+10*ISO*(1-AF); % [mM^-2 1/ms]
kom = 0.06; % [1/ms]
kiCa = 0.5; % [1/mM/ms]
kim = 0.005; % [1/ms]
ec50SR = 0.45; % [mM]
ks = 25; % [1/ms]

% Time discretization

% tspan = [0;1.54e5];
tspan = [0;100];

% tspan = [0;5e5];

%== Initial conditions

RyRr0 = 8.884332e-1;
RyRo0 = 8.156628e-7;
RyRi0 = 1.024274e-7;

y0 = [RyRr0 RyRo0 RyRi0];

% load('sternryr'); % load output of previous simulation saved as yfinal.mat
% y0 = yfinal;

% ODE15S
options = odeset('RelTol',1e-5,'MaxStep',1,'Stats','on'); 
[t,y] = ode15s(@(t,y) f(t,y,cajclamp,casrclamp,koCa,kom,kiCa,kim,ec50SR,ks), tspan, y0, options);

% yfinal = y(end,:);
% output = yfinal;
% save 'sternryr' 'yfinal'

figure(1)
plot(t,y(:,1),'b','LineWidth',5)
hold on
plot(t,y(:,2),'g','LineWidth',5)
plot(t,y(:,3),'m','LineWidth',5)
legend('R','O','I')
ylabel('\fontsize{45} Probability')
xlabel('\fontsize{45} Time(ms)')
% xlim([27640 3.3*10^4])
set(gca,'fontsize',45)

% save sternryrvar.mat y t;  % Save all vectors in y matrix. 

% Function
function ydot = f(t,y,cajclamp,casrclamp,koCa,kom,kiCa,kim,ec50SR,ks)
ydot = zeros(size(y));

%%

% RyR (Shannon/Grandi 2011)

MaxSR = 15; 
MinSR = 1;
kCaSR = MaxSR-((MaxSR-MinSR)/(1+(ec50SR/casrclamp)^2.5));
koSRCa = koCa/kCaSR;%
kiSRCa = kiCa*kCaSR;
RI = 1-y(1)-y(2)-y(3);
ydot(1) = (kim*RI-kiSRCa*cajclamp*y(1))-(koSRCa*(cajclamp^2)*y(1)-kom*y(2));   % R
ydot(2) = (koSRCa*(cajclamp^2)*y(1)-kom*y(2))-(kiSRCa*cajclamp*y(2)-kim*y(3));% O
ydot(3) = (kiSRCa*cajclamp*y(2)-kim*y(3))-(kom*y(3)-koSRCa*(cajclamp^2)*RI);   % I

I_rel = ks*y(2)*(casrclamp-cajclamp);          % [mM/ms]