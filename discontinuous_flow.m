%% The effect of disconnected nanopores with the capillary network. The base case presented assumes that the nanopores feeding the capillary network is in the range of .5nm to 2.5 nm. The nanopores connectd to the capillary network act as independent bottles and do not provide continuous gas supply till infinity. As a result the flow rate in the capilalry network fed by the nanopores will decrease with time. The aim of the code is to understand the time until which the high flow rate is sustained. This will highlight the importance of nanopore connectivity to support the long term production in shales


% this code is for radius 10 nm of pores.
% and another simulation shuld be for 2 nm
% 
clear all
close all
clc

%% Define the Figure size, color of the flow rates amd the colormap
f=figure('Renderer', 'painters', 'Position', [200 200 380 280]);
te=13;
col=["red","black","blue","green"]
c = colormap(hsv(5));

%% Define the time step, total time in seconds, Capillary radius and length, nanopore radius and length.

dt=0.00001;  % The time step
time=1:1:400;% The time is in seconds 
a=zeros(1,5); 
R=4e-6; %Radius of the capillary
L=4000e-6; % Lengh of the capillary
l=60e-6; % Length of the nanopore 
nr=[0.5e-9;10e-9;1.5e-9;2e-9;2.5e-9]; % radius of the nanopores

%% Define Inlet Flow Rate and gradient in capillary(Qin,grad) and Pressure 
% in the capillary and nanopores (CP, NP). Consta calculates the number of 
% nanopores that will populate the outer surface of the capillary

unit=(2*nr)/(0.4);
consta=(0.4*3.14*nr.^3)/(R^3*l);
CP=6894.76*[1000;2000;3000;4000;5000];
NP=CP+4;
Qin=[5.02e-16;7.53e-16;1e-15;1.25e-15;1.5e-15];
grad=-Qin/((3.14*R^4)/(8*2e-5));

%% Parameters for micro-capillary and nano-pores. c1,c2, l1, l2 are the 
% constants defined and derived in the analytical solution in the paper. 
% y is the x axis setp size. Qf is the final flow rate from the other end of the capillary
y=unit(2):unit(2):0.004;   
ggg=zeros(1,length(y));
P=6894.76*5000;
Pp=(P+4)*ones(1,length(y));
dx=y(1);
c1=consta(2);
c2=c1/dx;
l1=(c1+sqrt(c1^2+4*c2))/2;
l2=(c1-sqrt(c1^2+4*c2))/2;
Q=zeros(1,length(y));
Qf=zeros(length(time),1);
mf=zeros(length(time),1);

%% flow in nano-pores. Qout is the flow from the other end of the capillary.
Q_out=zeros(length(y),length(time));
mass_out=zeros(length(y),length(time));
volume=3.14*(nr.^2)*l;
Q_initial=volume(2)*ones(length(y),length(time));
P_all=zeros(length(y),length(time));
P_all_diff=zeros(length(y),length(time));

%% density, mass and moles in nano-pores
density=0.158*ones(1,length(y));
mass=(volume.*(16.04/101.15));
initial_mass=mass(2)*ones(1,length(y));
left_mass=zeros(1,length(y));
moles=(volume.*(1/101.15));
initial_moles=moles(2)*ones(1,length(y));
left_moles=zeros(1,length(y));
ratio=zeros(1,length(time));

%% to get the no of nano-pores around the circumference of the capillary
mul=zeros(length(nr),1);
for i=1:length(nr)
    mul(i)=(0.4*2*3.14*R)/(2*nr(i));
end
mult=(3.14*nr.^4)/(8*2e-5*l);

%% simulation starts
for ttt=2:5  
    P=6894.76*5000;
    Pp=(P+4)*ones(1,length(y));
    Q=zeros(1,length(y));
    Qf=zeros(length(time),1);
    mf=zeros(length(time),1);
    % flow in nano-pores
    Q_out=zeros(length(y),length(time));
    mass_out=zeros(length(y),length(time));
    volume=3.14*(nr.^2)*l;
    Q_initial=volume(2)*ones(length(y),length(time));
    P_all=zeros(length(y),length(time));
    P_all_diff=zeros(length(y),length(time));
    % density, mass and moles in nano-pores
    density=0.158*ones(1,length(y));
    mass=(volume.*(16.04/101.15));
    initial_mass=mass(2)*ones(1,length(y));
    left_mass=zeros(1,length(y));
    moles=(volume.*(1/101.15));
    initial_moles=moles(2)*ones(1,length(y));
    left_moles=zeros(1,length(y));

    step_wise_Pp=zeros(length(y),length(time));
    step_wise_gg=zeros(length(y),length(time));  
    for m=1:length(time)
    
        for b=1:length(y)
            ggg(1,b)=((P-Pp(1,b))/(l1-l2) )*(l1*exp(l2*y(b))-l2*exp(l1*y(b)));
            ggg(1,b)=ggg(1,b)+ (-grad(ttt)/(l1-l2))*(exp(l2*y(b))-exp(l1*y(b))) +Pp(1,b);
        end
        di=zeros(1,length(y));
        
        for i = 1:length(y)
            di(1,i)=ggg(1,i)-P;
        end
        Pressure_diff=zeros(1,length(y));
        for j=1:length(y)
            Pressure_diff(1,j)=Pp(1,j)-ggg(1,j);
        end
        % to get the flow from each nano-pore
        for k=1:length(y)
            Q(1,k)=mult(2)*Pressure_diff(1,k);
        end        
        % gas coming out from the nano-pore in time dt
        for o=1:length(Q)
            Q_out(o,m)=Q(1,o)*dt;
            mass_out(o,m)=Q(1,o)*dt*density(1,o);
        end        
        step_wise_Pp(1:length(y),m)=Pp(1,:);
        step_wise_gg(1:length(y),m)=ggg(1,:);           
        % to get the gas left in the nanopore after time dt
        for h=1:length(y)
            Q_initial(h,m:length(time))=Q_initial(h,m:length(time))-Q_out(h,m);
        end
        % to get Qout ater time dt
        Qf(m,1)=sum(Q_out(:,m));
        mf(m,1)=sum(mass_out(:,m));
        
        % New Pressure, volume, density in the nano-pore
        for d=1:length(y)
            left_mass(1,d)=initial_mass(1,d)-Q_out(d,m)*(0.158);
        end
        
        for e=1:length(y)
            left_moles(1,e)=initial_moles(1,e)-Q_out(e,m)*(1/101.15);
        end
        
        % Final Pressure in the nano-pores.
        for f=1:length(y)
            Pp(1,f)=Pp(1,f)*(left_moles(1,f)/initial_moles(1,f)); 
        end
       
        tot=mf(m,1)*mul(2);
        f=tot+(Qin(ttt)*dt*0.158);
        ratio(1,m)=f/(Qin(ttt)*dt*0.158);
        
        initial_moles=left_moles;
        initial_mass=left_mass;
    end        
plot(time/100,ratio,'linewidth',1.5,'color',col(ttt-1));
hold on
end

%% plotting the graph

grid on
box on
ax = gca; % current axes
ax.FontSize = 14;
ax.TickDir = 'in';
ax.GridAlpha = 0.13;
ax.FontWeight = 'normal';
xlabel('$$\mathrm{Time,~\mathrm{ms}}$$','interpreter','latex')
ylabel('$$\mathrm{Ratio~of~flow~rates,~Q_{out}/Q_{in}}$$','interpreter','latex')
s = ...
    {'$$Q_{in} =  753~{\mu}\mathrm{m}^3\mathrm{/s}$$';
     '$$Q_{in} = 1000~{\mu}\mathrm{m}^3\mathrm{/s}$$';
     '$$Q_{in} = 1250~{\mu}\mathrm{m}^3\mathrm{/s}$$';
     '$$Q_{in} = 1500~{\mu}\mathrm{m}^3\mathrm{/s}$$'};
legend(s,'Location','northeast','interpreter','latex')
set(gca,'LineWidth',0.2,'TickLength',[0.007 0.007]);
print('-depsc2','-r400','effect_of_time.eps');
