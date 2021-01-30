%plotted all particles on the same graph 
%not enough time to seperate due to other assignmen due
close all;
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%electron spec
 global C

    addpath ../geom2d/geom2d

    C.hb = 1.054571596e-34;             % Dirac constant
    C.h = C.hb * 2 * pi;                    % Planck constant
    C.m_0 = 9.10938215e-31;             % electron mass
    C.kb = 1.3806504e-23;               % Boltzmann constant
    C.eps_0 = 8.854187817e-12;          % vacuum permittivity
    C.mu_0 = 1.2566370614e-6;           % vacuum permeability
    C.c = 299792458;                    % speed of light
    C.g = 9.80665; %metres (32.1740 ft) per s²

nTime = 100; %number of steps
dt = 1;
%EFx = 5 *10^-20; %electric field in x
force = 1;
%deltaT = 1E-6;
%T = 0;
%mass = 1E-31;
nPoints=100;
nSims = 100;
n=0;
xInitial=0
xPosition = zeros(1,2);
Velocity = zeros(1,2);
%Ps = 0.05;
Ps = 1-exp(-dt/20)
driftV = zeros(1,2);
t=0;

%calutalted acceleration due to force
    for i = 2: nSims
         for j = 1:2
         %drift_velocity = zeros(num_of_eletr, points_num);
         %Avg = mean()
         
         prob = rand(1)*1;

         %to reinitialize if in range
         if(prob<=Ps)
             Velocity(i,j) = 0;%starting v
         else 
             Velocity(i,j) = Velocity(i-1,j) + ((C.q_0*force)/C.m_0)*dt;%calculating v
         end
         
         driftV(i) = mean(Velocity(i));    %drift calculation
         overallDrift(j) = mean(driftV);
         xPosition(i,j) = xPosition(i-1,j) + Velocity(i-1,j)*dt;      %position array
         t(i)= t(i-1)+1;    %time array
         
         title(sprintf("V vs time. Average =%0.5f",driftV(i)));
         subplot(3,1,1)
         plot(t,Velocity,'b')
         plot(t,driftV(i),'ro')
         drawnow
         hold on
         
         title("V vs position");
         subplot(3,1,2)
         plot(xPosition,Velocity,'r')
         plot(xPosition,driftV,'ro')
         drawnow
         hold on
         
         title("Position vs time");
         subplot(3,1,3)
         plot(t,xPosition,'g')
         drawnow
         hold on
           
         end
    end

