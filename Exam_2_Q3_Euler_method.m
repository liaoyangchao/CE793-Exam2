function Yangchao_Euler_method()
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Exam 2
%%%% Student Name: Yangchao Liao
%%%% Student ID.: 1299252
%%%% Department: Civil & Environmental Eng.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
close all;
clc;

%% Initial and boundary conditions
Nx = 10;
Lx = 15;

dx = Lx/(Nx-1);
x = 0:dx:Lx;
alpha = 1;

T = zeros(1,Nx);

T(1) = 0;     % T(0,t) = 0
T(Nx) = Lx^2*exp(-Lx);   % T(Lx,t) = Tsteady(Lx)

dt = 0.001;

%% Explicit Euler solution

T_Euler = zeros(1,Nx);
T_steady_check = zeros(1,Nx);

% Time loop 
for iteration = 1:100000
  iteration;
for i=2:Nx-1
    T_Euler(i) = T_Euler(i) + dt * ( alpha * (T_Euler(i+1) - 2*T_Euler(i) + T_Euler(i-1))/(dx^2)...
        - (x(i)^2 - 4*x(i) + 2)*exp(-x(i)) );
    
    % Steady state checking
    T_steady_check(i) =  alpha * (T_Euler(i+1) - 2*T_Euler(i) + T_Euler(i-1))/(dx^2)...
        - (x(i)^2 - 4*x(i) + 2)*exp(-x(i));
    
end
     T_Euler(floor(Nx/2));
     T_steady_check(floor(Nx/2));
     max(max(abs(T_steady_check)))
end

figure(1)
P = plot(x,T_Euler,'-','LineWidth',2.5,'MarkerFaceColor','w','MarkerSize',10);hold on


%% Initial and boundary conditions
Nx = 20;
Lx = 15;

dx = Lx/(Nx-1);
x = 0:dx:Lx;
alpha = 1;

T = zeros(1,Nx);

T(1) = 0;     % T(0,t) = 0
T(Nx) = Lx^2*exp(-Lx);   % T(Lx,t) = Tsteady(Lx)

dt = 0.001;

%% Explicit Euler solution

T_Euler = zeros(1,Nx);
T_steady_check = zeros(1,Nx);

% Time loop 
for iteration = 1:100000
  iteration;
for i=2:Nx-1
    T_Euler(i) = T_Euler(i) + dt * ( alpha * (T_Euler(i+1) - 2*T_Euler(i) + T_Euler(i-1))/(dx^2)...
        - (x(i)^2 - 4*x(i) + 2)*exp(-x(i)) );
    
    % Steady state checking
    T_steady_check(i) =  alpha * (T_Euler(i+1) - 2*T_Euler(i) + T_Euler(i-1))/(dx^2)...
        - (x(i)^2 - 4*x(i) + 2)*exp(-x(i));
    
end
     T_Euler(floor(Nx/2));
     T_steady_check(floor(Nx/2));
     max(max(abs(T_steady_check)))
end

figure(1)
P1 = plot(x,T_Euler,'-','LineWidth',2.5,'MarkerFaceColor','w','MarkerSize',10);hold on

%% Exact solution

Nx = 100;
Lx = 15;

dx = Lx/(Nx-1);
x = 0:dx:Lx;
alpha = 1;

T_steady = x.^2 .* exp(-x);

%% Plotting T(x)

P2 = plot(x,T_steady,'-','LineWidth',2.5,'MarkerFaceColor','w','MarkerSize',10);hold on

lgd=legend('Explicit Euler method: N_x = 10','Explicit Euler method: N_x = 20','Exact solution');
set(lgd,'Location','NorthEast','Orientation','vertical');set(lgd,'Box','off'); %vertical

xlabel('x','FontName','Arial','FontSize',25)
ylabel('T','FontName','Arial','FontSize',25)

set(gca,'linewidth',1.5,'FontName','Arial','FontSize',25);
set(gcf,'Color','w','Units','inches','position',[0,0,8,6]);
AxesH = gca;InSet = get(AxesH, 'TightInset');
set(AxesH, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
hold off;
end
