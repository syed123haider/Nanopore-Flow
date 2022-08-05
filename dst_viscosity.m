%%%%%%%%%%%%%%%%%%%%%%%%%%%%calculates density profile%%%%%%%%%%%%%%%%%%%
clear all;clc;
close all;
%Tf=176;
T=373;
%T= (Tf + 459.67)*(5/9);
P=6894.76*5000;
Tc=190;
Pc=4.600183872000000e+06;
w=0.011;
MW=16/1000;
Liquido=0;
R = 8.314; % gas constant [=] J/(mol K)
 
% Reduced variables
Tr = T/Tc ;
Pr = P/Pc ;
 
% Parameters of the EOS for a pure component
m = 0.37464 + 1.54226*w - 0.26992*w^2;
alfa = (1 + m*(1 - sqrt(Tr)))^2;
a = 0.45724*(R*Tc)^2/Pc*alfa;
b = 0.0778*R*Tc/Pc;
A = a*P/(R*T)^2;
B = b*P/(R*T);
 
Z = roots([1 -(1-B) (A-3*B^2-2*B) -(A*B-B^2-B^3)]);
 
ZR = [];
for i = 1:3
   if isreal(Z(i))
    ZR = [ZR Z(i)];   
   end
end

if Liquido == 0
    Z = min(ZR);   
else
    Z = max(ZR);
end
 
% Fugacity coefficient
fhi = exp(Z - 1 - log(Z-B) - A/(2*B*sqrt(2))*log((Z+(1+sqrt(2))*B)/(Z+(1-sqrt(2))*B)));
if isreal(fhi)
    density=P*MW/(Z*R*T);
    result = [Z fhi density];
else
    'No real solution for "fhi" is available in this phase'
    result=['N/A' 'N/A' 'N/A'];
end
 
den=(result(1,3))/1000;
den=0.226;
base=den;


%%%%%%%%%%%%%%%% Bulk fugacity calculation %%%%%%%%%%%%%%%%%%%%%%
 
den=(den*1000000)/16;
ln_fb= (b*den)/(1-b*den) - ((a*den)/(R*T))*(1/(1+2*b*den-b^2*den^2));
ln_fb= ln_fb -  log((1-b*den)/(R*T*den)) - (a/(2*sqrt(2)*R*T*b))*log((1+(1+sqrt(2))*den*b)/(1+(1-sqrt(2))*den*b));
 
fb=exp(ln_fb);
 
%%%%%%%%%%%%%%%% si calculation %%%%%%%%%%%%%%
Na=6.023e23;
sigma_f=0.3758e-9;
sigma_s=0.34e-9;
sigma_fs=(sigma_f+sigma_s)/2;
rho_atom=38.2e18;
e_f=148.6;
e_s=28;
kb=1.38064852e-23;
e_fs=sqrt(e_s*e_f)*kb;
%%%% definition of the slit


L=[0.4e-9;0.5e-9;0.6e-9;0.7e-9;0.8e-9;0.9e-9;1e-9;2e-9];
% L=[0.4e-9;0.5e-9;0.6e-9;0.7e-9;0.8e-9;0.9e-9;1e-9;2e-9;2e-9;4e-9;6e-9;8e-9;10e-9];
avg_den=zeros(1,length(L));
density_data=zeros(2/0.005,length(L));
distance_data=zeros(2/0.005,length(L));
ttt=0;
f=figure('Renderer', 'painters', 'Position', [200 200 380 280]);
te=13;
for g=1:length(L)
    l=L(g);
    l
    new=0.18e-9:0.005e-9:l-0.18e-9;
    distance_data(1:length(new),g)=new;
    final_d=zeros(1,length(new));
    
        for j=1:length(new)
            z=new(j);
            ttt=ttt+1;
            z_d1=z+(sigma_f)/2;
            z_d2=L(g)-z + sigma_f/2;
            
            sum1=0;
            sum2=0;
                for i=1:4
                     sum1=sum1+(sigma_fs^4)/((z_d1+(i-1)*sigma_s)^4);
                end
 
                for i=1:4
                     sum2=sum2+(sigma_fs^4)/((z_d2+(i-1)*sigma_s)^4);
                end
                
                
        si_1=4*3.14*rho_atom*e_fs*sigma_fs^2*((sigma_fs^10)/(5*z_d1^10) - 0.5*sum1);
        si_2=4*3.14*rho_atom*e_fs*sigma_fs^2*((sigma_fs^10)/(5*z_d2^10) - 0.5*sum2);
        mu_fs=Na*(si_1+si_2);
        f_ff=fb*exp(-(mu_fs)/(R*T));
        
    if (L(g)/sigma_f >=3)
        
        if z/sigma_f<0.5
        ad=0.5 +5/6 -(1/3)*(((L(g)/sigma_f-0.5)-0.5)^-3); 
        ad=ad*a*(3/8);
        ex=11
        elseif z/sigma_f>=0.5 & z/sigma_f<=1.5
        ad=z/sigma_f +5/6 -(1/3)*((((L(g)-z)/sigma_f)-0.5)^-3); 
        ad=ad*a*(3/8);
        ex=12
        elseif z/sigma_f>=1.5 & z/sigma_f <= L(g)/sigma_f-1.5
        ad=8/3-1/3*(((z/sigma_f)-0.5)^-3) -(1/3)*(((((L(g)-z)/sigma_f)) -0.5)^-3);
        ad=ad*a*(3/8);
        ex=13
        elseif z/sigma_f >= L(g)/sigma_f-1.5 & z/sigma_f <= L(g)/sigma_f -0.5
        ad= ((L(g)-z)/sigma_f) +5/6 - (1/3)*(((z/sigma_f)-0.5)^-3);
        ad=ad*(3/8)*a;
        ex=14
        elseif z/sigma_f>L(g)/sigma_f-0.5
        ad=0.5 +5/6 -(1/3)*(((L(g)/sigma_f-0.5)-0.5)^-3); 
        ad=ad*a*(3/8);
        ex=15
        end
     end
 
 
    if (L(g)/sigma_f >2 & L(g)/sigma_f<3)
        
        if z/sigma_f<0.5
        ad=0.5+5/6 -(1/3)*(((L(g)/sigma_f-0.5)-0.5)^-3); 
        ad=ad*a*(3/8);
        ex=20
        elseif z/sigma_f>=0.5 & z/sigma_f<= L(g)/sigma_f - 1.5
        ad=z/sigma_f +5/6 -(1/3)*((((L(g)-z)/sigma_f)-0.5)^-3); 
        ad=ad*a*(3/8);
        ex=21
        elseif z/sigma_f >= L(g)/sigma_f-1.5 & z/sigma_f <= 1.5
        ad=(L(g)/sigma_f)-1;
        ad=ad*(3/8)*a;
        ex=22
        elseif z/sigma_f >= 1.5 & z/sigma_f <= L(g)/sigma_f -0.5
        ad= ((L(g)-z)/sigma_f) +5/6 - (1/3)*(((z/sigma_f)-0.5)^-3);
        ad=ad*(3/8)*a;
        ex=23
        elseif z/sigma_f > L(g)/sigma_f-0.5
        ad=0.5+5/6 -(1/3)*(((L(g)/sigma_f-0.5)-0.5)^-3); 
        ad=ad*a*(3/8);
        ex=24
        end
    end
 
    if (L(g)/sigma_f >=1.5 & L(g)/sigma_f <=2)
        ad=a*(3/8)*(L(g)/sigma_f-1);
        ex=31
    end
 
    if (L(g)/sigma_f >1 & L(g)/sigma_f <= 1.5)
        ad=a*(3/16)
        ex=41
    end
    
    d=0:1:40000;
    val=zeros(1,length(d));
        for i=1:length(d)
            val(i) = log(f_ff)-(b*d(i))/(1-b*d(i))+((ad*d(i))/(R*T))*(1/(1+2*b*d(i)-b^2*d(i)^2))+log((1-b*d(i))/(R*T*d(i))) + (ad/(2*sqrt(2)*R*T*b))*log((1+(1+sqrt(2))*d(i)*b)/(1+(1-sqrt(2))*d(i)*b));
        end
        
        brea=min(val);
        count=0;
            for i=1:length(val)
                check=val(i);
                if (check == brea)
                    break
                end
                count=count+1;
            end
            final_d(j)=d(count);
            min(val)
        end
       data=final_d;
       data=data.*(16/1000000);
       density_data(1:length(data),g)=final_d;
       avg_den(g)=sum(final_d)*((0.005e-9*16)/(1000000*(L(g)-2*0.18e-9)));
end

avg_den=avg_den';
avg_den=avg_den./base;
base=ones(1,length(new))*base;

%%%%%%%%%%%%%%%%%%%%%%viscosity calculation%%%%%%%%%%%%%%%%%%
dimm=size(density_data)
c = colormap(hsv(dimm(2)));
a_Tc=(0.45724*R^2*Tc^2)/(Pc);
m=0.37464+1.54226*w-0.26997*w^2;
alpha=(1+m*(1-Tr^0.5))^2;
a=a_Tc*alpha;
bdash=0.07780*R*Tc/Pc;

for i=1:dimm(2)
    count=0;
    for j=1:length(distance_data(:,i))
           if (distance_data(j,i) == 0)
                break
           end
         count=count+1;
     end
    X=distance_data(1:count,i);
    Y=density_data(1:count,i);
    V=1./Y;
    visc=zeros(length(V),1);
    for o=1:length(V)
        brx= (R/(V(o)-bdash)) + ((a_Tc*(m+m^2*(1-Tr^0.5)))/(V(o)*(V(o)+bdash)+bdash*(V(o)-bdash)))*(1/(T*Tc)^0.5);
        brx=brx*(V(o)/R);
        br=(bdash+ a_Tc*(m+m^2*(1-Tr^0.5))/(R*sqrt(T*Tc)))*(1/V(o));
        b=(bdash+ a_Tc*(m+m^2*(1-Tr^0.5))/(R*sqrt(T*Tc)));
        x=brx/br;
        ccc=br*((1/brx)+0.8+0.7614*(brx)^2);
        sigma=0.3758e-9;
        k=1.38064852e-23;
        m1=16/(1000*6.023e23);
        n0=(1.016*(5/16)*(3.14*m1*k*T)^0.5)*(1/(3.14*sigma^2));
        visc(o)=n0*ccc;
    end
    plot(X*1000000000,visc*1000,'linewidth',1.5,'color',c(i,:))
    hold on 
end
bulk_den=(base*1000000)/16;
V=1/bulk_den(1);
brx= (R/(V-bdash)) + ((a_Tc*(m+m^2*(1-Tr^0.5)))/(V*(V+bdash)+bdash*(V-bdash)))*(1/(T*Tc)^0.5);
brx=brx*(V/R);
br=(bdash+ a_Tc*(m+m^2*(1-Tr^0.5))/(R*sqrt(T*Tc)))*(1/V);
b=(bdash+ a_Tc*(m+m^2*(1-Tr^0.5))/(R*sqrt(T*Tc)));
x=brx/br;
ccc=br*((1/brx)+0.8+0.7614*(brx)^2);
sigma=0.3758e-9;
k=1.38064852e-23;
m1=16/(1000*6.023e23);
n0=(1.016*(5/16)*(3.14*m1*k*T)^0.5)*(1/(3.14*sigma^2));
visc=n0*ccc;
visc=visc*ones(length(X),1);

plot(X*1000000000,visc*1000,'linewidth',1.5);

grid on
box on
ax = gca; % current axes
ax.FontSize = 14;
ax.TickDir = 'in';
ax.FontWeight = 'normal';
ax.GridAlpha = 0.13;
%xlim([0 1])
%ylim([0 0.55])
xlabel('$$\mathrm{Pore~Diameter,nm}$$','interpreter','latex')
ylabel('$$\mathrm{Methane~viscosity,(cP)}$$','interpreter','latex')
s = ...
    {
     '$$0.4~\mathrm{nm}$$'
     '$$0.5~\mathrm{nm}$$';
     '$$0.6~\mathrm{nm}$$';
     '$$0.7~\mathrm{nm}$$';
     '$$0.8~\mathrm{nm}$$';
     '$$0.9~\mathrm{nm}$$';
     '$$1~\mathrm{nm}$$';
     '$$2~\mathrm{nm}$$';
     '$$\mathrm{Bulk~viscosity}$$'};
legend(s,'Location','northeast','interpreter','latex')
set(gca,'LineWidth',0.2,'TickLength',[0.007 0.007]);
print('-depsc2','-r400','dst_viscosity.eps');