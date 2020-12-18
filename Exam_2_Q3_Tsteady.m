function Yangchao_Tsteady()
  
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
Nx = 100;
Lx = 15;

dx = Lx/(Nx-1);
x = 0:dx:Lx;
alpha = 1;

T_steady = x.^2 .* exp(-x);

%% Plotting T(x)
figure(1)
P = plot(x,T_steady,'-','LineWidth',2.5,'MarkerFaceColor','w','MarkerSize',10);hold on


xlabel('x','FontName','Arial','FontSize',25)
ylabel('T_{steady}','FontName','Arial','FontSize',25)

set(gca,'linewidth',1.5,'FontName','Arial','FontSize',25);
set(gcf,'Color','w','Units','inches','position',[0,0,8,6]);
AxesH = gca;InSet = get(AxesH, 'TightInset');
set(AxesH, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
hold off;
end
