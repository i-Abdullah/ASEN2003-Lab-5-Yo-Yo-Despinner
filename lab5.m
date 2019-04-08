%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Lab 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Housekeeping and Constants
clc; clear all; close all;

omega0 = 130; % RPM
I = 0.0063; %+/- 0.0001 kg*m^2
R = 0.0762; %m
m = .054*2; %kg

%% Theoretical Tangential
omega0 = omega0*2*pi/60; %rad/s
C = I/(m*R^2)+1;
t = linspace(0,2,100);
a = (C+((omega0^2)*(t.^2)));
b = (C-((omega0^2)*(t.^2)));
omega = omega0*(b./a);
alpha = (-4.*C.*(omega0^3).*t)./((C+(omega0^2).*t.^2).^2);
T = m*((omega.*R).^2)/R; 
subplot(3,1,1)
plot(t,omega)
ylabel('Omega (rad/s)')
title('Theoretical Kinematic and Force Values with Tangential Release');
subplot(3,1,2)
plot(t,alpha)
ylabel('Alpha (rad/s^2)')
subplot(3,1,3)
plot(t, T)
ylabel('Tension Force (N)')
xlabel('Time (s)')
Lt = R*sqrt(C); %m tangential
Lr = R*sqrt((I/(m*R^2)+1)-1); %Given
Ltest = sqrt(I/m+R^2)-R; %Derived


%% Experimental omega
data1 = load('YoYo_despin');
data2 = load('YoYo_noweight');
timeDe = data1(:,1);
omegaDe = data1(:,2);




omegaDe = omegaDe.*(0.104719755);%rpm to rad/s
timeNo = data2(:,1);
omegaNo = data2(:,2);
omegaNo = omegaNo.*(0.104719755);%rpm to rad/s


% calculate avg omega0 that the experiment started with so the model starts
% with too!

omega0_avg_exp = mean(omegaDe(1:18));


%Fixing errors in time data
for i = 120:267
    timeDe(i) = timeDe(i)+1180;
end
for i = 268:1001
    timeDe(i) = timeDe(i)+1180+1480;    
end
for i = 126:275
    timeNo(i) = timeNo(i)+1240;
end
for i = 276:3922
    timeNo(i) = timeNo(i)+1240+1500;
end

% zero the matrices:
timeDe(1) = [];
omegaDe(1) = [];

% eyeball where we started the mechanism, then zero from there:
timeDe(1:19) = [];
omegaDe(1:19) = [];

% zero time:

timeDe = timeDe - timeDe(1);

figure(2)
plot(timeNo./1000, omegaNo,'LineWidth',1.5)
legend('No Masses Used')
xlabel('Time (s)')
ylabel('Omega (rad/s)')
title('Angular Velocity Without activating the YOYO mechanism')
grid minor
ylim([ -2 15])
xlim([ -1 42])

% estimate friction moment:

AlphaNo = diff(omegaNo(2:end));

MomentFriciton = (sum(AlphaNo(2:end))/length(AlphaNo(2:end)))/(sum(diff(timeDe(2:end)))/length(timeDe(2:end)))*I ;

MomentFriciton = (-I*omegaNo(2)) / (timeNo(end)./1000) ;


%% Experimental alpha

alphaDe = diff(omegaDe)./diff(timeDe);
alphaNo = diff(omegaNo)./diff(timeNo);

%MomentFriciton = mean(alphaNo(2:end))*I ;


omega = @(t) omegaDe(1)*(((C-((omega0^2)*(t.^2))))./((C+((omega0^2)*(t.^2)))));

TimePlug = (timeNo./1000);


% compute when tangential realse ends:


Wtang = omegaDe(1)*((C-((omegaDe(1)^2)*(Lr.^2)))/(C+((omegaDe(1)^2)*(Lr.^2))));


timeDe = timeDe(1:end);
timeNo = timeNo(1:end);
figure(3)
plot(timeDe(1:end)./1000,omegaDe(1:end),'-','LineWidth',1.5)
hold on
plot(TimePlug(2:end),omega(TimePlug(2:end)),'LineWidth',1.5)
plot(0.185,3.5,'ok','MarkerSize',8)
legend('Despin Masses Used (radial)','No Masses Used (theoretical, tangential)','Point when tangential release ends')
ylabel('\omega (rad/s^2)')
xlabel('Time (s)')
title('Angular Velocity theoretical vs experimental')
xlim([-0.01 1.5])
grid minor


