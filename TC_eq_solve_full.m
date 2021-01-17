%Script to find how the pressures at the front of the detector vary when
%the signal changes

%Initialise the variables
clear p1 p2 Q S_p S_i S_s V_1 V_2

reset(symengine)

syms p1(t) p2(t) Q S_p S_i S_s V_1 V_2

%Setup the differential equations
ode1 = V_1*diff(p1)== Q - S_p*(p1-p2)-S_s*p1;
ode2 = V_2*diff(p2) == S_p*(p1-p2) - S_i*p2;

odes=[ode1;ode2];

%Dp1=diff(p1);
%Dp2=diff(p2);

%Give the initial conditions
cond1 = p1(0) == 0;
cond2 = p2(0) == 0;

conds=[cond1;cond2];

%Solve the equations
S = dsolve(odes,conds);

p1_sol=S.p1;
p2_sol=S.p2;

%substitute for values
d=16e-3; %Diameter of CF16 in m
L=98e-3+45e-3;%90e-3+120e-3+25e-3; %Length of CF16 from the sample chamber to the detector. +25e-3 due to T piece volume sticking up.
r_ap=0.05; %Radius of detector aperture in cm 


V_1=pi*(d/2)^2*L*1000; %Volume of the input pipe to the detector
V_2=(3706.9e-6+pi*(2.5e-3^2+8.5e-3^2))*35.8e-3*1000;
%V_2=pi*(33.5e-3)^2*35.8e-3*1000; %Lower bound on the volume inside the front of the detector
V_3=pi*(28e-3/2)^2*68e-3*1000;
V_2=V_2+V_3;

%Quantites to calculate average speed
k_B=1.38064852e-23;
T=293;
m=4*1.6726219e-27;

c_bar=sqrt((8*k_B*T)/(pi*m)); %Average speed of atoms
c_bar_cm=c_bar*100; %Average speed in cm

C_pipe_pre=c_bar_cm*pi/12*1e-3; %Prefactor for conductance of pipe
C_ap_pre=c_bar/4*0.1; %Prefactor for conductance of aperture

S_p=C_pipe_pre*(0.5)^3/(3.884); %Pumping speed of the small pipe into the detector
S_i=0.67;%Estimate of the pumping speed down the ioniser
%S_i=C_pipe_pre*(1.22)^3/(66.1);
%C_pipe_pre*(2.8)^3/(7.6);

S_s=C_ap_pre*(pi*r_ap^2); 

S_out=(1/S_p+1/S_i)^-1;%Pumping speed to the pump from the front
Stot=S_out+S_s; %Total pumping speed in front of ioniser

Q=1;

p1_lim=Q/Stot;
%p2_lim=Q/S_i;
p2_lim=Q/(S_i+S_s*(1+S_i/S_p));

%Plot results
p1_sol_sub=subs(p1_sol);
p2_sol_sub=subs(p2_sol);

p1_h=matlabFunction(p1_sol_sub);
p2_h=matlabFunction(p2_sol_sub);

t_plot=0:1e-3:2;
figure;
plot(t_plot,p1_h(t_plot),'b')
hold on
plot([min(t_plot), max(t_plot)],[p1_lim,p1_lim],'b--','Linewidth',1)
plot(t_plot,p2_h(t_plot),'r')
plot([min(t_plot), max(t_plot)],[p2_lim,p2_lim],'r--','Linewidth',1)
legend('p_1 - Input pipe pressure','Expected p_1 limit','p_2 - Pressure at front of the ioniser','Expected p_2 limit','Location','SouthEast')

xlabel('Time/s')
ylabel('Pressure/a.u.')
set(gca,'FontSize',12)
set(gca,'Linewidth',1)

%print -depsc2 C:\Users\Matt\Dropbox\Thesis\Detector\Figures\TC_eqs.eps

P1_out=p1_h(t_plot);
P2_out=p2_h(t_plot);

Det_eff=S_out/(S_out+S_s);
%Random variables, ignore.
%a=S_i^2*V_1^2 + 2*S_i*S_p*V_1^2 - 2*S_i*S_p*V_1*V_2 + S_p^2*V_1^2 + 2*S_p^2*V_1*V_2 + S_p^2*V_2^2;
%b2=((S_i*V_1 + S_p*V_1 + S_p*V_2 - (a)^(1/2)))/(2*V_1*V_2);
% % Fit: 'untitled fit 1'.
% [xData, yData] = prepareCurveData( t_plot, P2_out );
% 
% % Set up fittype and options.
% ft = fittype( 'a*exp(-(x/b))+c', 'independent', 'x', 'dependent', 'y' );
% excludedPoints = excludedata( xData, yData, 'Indices', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46] );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.StartPoint = [0.278498218867048 0.546881519204984 0.957506835434298];
% opts.Exclude = excludedPoints;
% 
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
% 
% error_mat=confint(fitresult);
% error=(error_mat(2,2)-error_mat(1,2))/2;


%Solve the diffusion equation numerically

N_t=1000000;
N_x=100;

t_end=0.5;
C_sol=0.67/1000; %Convert to m^3
L=0.661;%0.8;
A=pi*(12.2e-3/2)^2;

%Diffusion contant
D=C_sol*L/A;

p=NaN*zeros(N_x,N_t);

delta_t=t_end/N_t;
delta_x=L/N_x;

t_sol=linspace(0,2,N_t);

%Initial conditions at t=0
p(:,1)=0;

%Initial conditions at x=0
p_interp=interp1(t_plot,P2_out,t_sol);
p(1,:)=p_interp;

%Initial conditions at x=L
p(end,:)=0;

%Iterate each time step and update the new pressures

for n=1:(N_t-1)
    for m=2:(N_x-1)
        p(m,n+1)=p(m,n)+D*delta_t/(delta_x^2)*(p(m+1,n)-2*p(m,n)+p(m-1,n));
    end
end


p_int=sum(p,1)*delta_x;

t_reduce=t_sol(1:100:end);
p_reduce=p_int(1:100:end);

%figure;plot(t_sol,p_int);


clip_time=0.1;

%% Fit an exponential to the obtained signal
[xData, yData] = prepareCurveData( t_reduce, p_reduce );

% Set up fittype and options.
ft = fittype( 'a*exp(-x/b)+c', 'independent', 'x', 'dependent', 'y' );
excludedPoints = excludedata( xData, yData, 'Domain', [clip_time 2] );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.679702676853675 0.655098003973841 0.162611735194631];
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

clip_ind=find(t_reduce<clip_time, 1, 'last' );

figure;plot(t_reduce(1:clip_ind),p_reduce(1:clip_ind),'r.')
hold on
plot(t_reduce(clip_ind+1:end),p_reduce(clip_ind+1:end),'k.')
plot(t_reduce,fitresult(t_reduce),'Color',[0,0.4470,0.7410],'LineWidth',2) %[0.8500,0.3250,0.0980]
ylim([0 0.35])
xlim([0 1])
set(gca,'FontSize',12,'LineWidth',1)
xlabel('Time/s')
ylabel('Ion current/a.u.')

error_mat=confint(fitresult);
error=(error_mat(2,2)-error_mat(1,2))/2;

% print -depsc2 C:\Users\Matt\Dropbox\Thesis\Detector\Figures\TC_sim.eps
% savefig(['C:\Users\Matt\Dropbox\Thesis\Detector\Figures\TC_sim.fig'])

