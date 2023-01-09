%% Gamma distributed DDE simulation

% Thi simulation produces the simulation plots for the simulation dynamics corresponding to the
% convergence plots of the the FCRK method  on the log-log scale. For
% the linear and non-linear tests, we must have j as an integer to produce
% a reference solution through the linear chain technique.

close all
clear all

% The time interval for simulation and fixed stepsize hIn.
TMin = 0;
TMax =  100;
hIn = 1e-2;

%% Linear test problem
% Model parameters
PA.Alpha =  0.8; 
PA.Beta =    -11/10;

% Delay distribution parameters
tau =  1.0; % 1.5;
j = 1;   

%Right handside of the distributed DDE where $d$ is the delayed argument sand phi is the history function
F = @(x,d) PA.Alpha*x+PA.Beta*d;
phi = @(t) 1 ;

[solDDE1] = ddef4(F, phi, tau, j, TMin, TMax, hIn);

j = 4;
[solDDE2] = ddef4(F, phi, tau, j, TMin, TMax, hIn);

j = 7;
[solDDE3] = ddef4(F, phi, tau, j, TMin, TMax, hIn);

Points = [TMin:hIn:TMax]';

%% Figures
Fig1 = figure(1);
g1 = plot(solDDE1.x,solDDE1.y(:,1),'LineWidth',1.75,'Color', [171,217,233]/255,'LineStyle','-');
hold on 
g2 = plot(solDDE2.x,solDDE2.y(:,1),'LineWidth',1.75,'Color', [118,42,131]/255,'LineStyle','--');
hold on 
g3 = plot(solDDE3.x,solDDE3.y(:,1),'LineWidth',1.75,'Color', [90,174,97]/255,'LineStyle',':');
hold on 
ylabel('y(t)','FontSize',15); % ,'Interpreter','latex','FontSize',15)
xlabel('Time','FontSize',15)
legend([g1 g2 g3],'j = 1','j = 4', 'j = 7', 'Location','SouthEast') 


%%  Nonlinear test problem
%Carrying capacity
% PA.K = 2;

%Right handside of the distributed DDE where $d$ is the delayed argument and phi is the history function

% F = @(x,d) x-x*d/PA.K; % +cos(d);
% phi = @(t) 1; % 1+sin(t);

% tau =  2.25; % 1.25;
% j = 3;  % 
% 
% [solDDE1] = ddef4(F, phi, tau, j, TMin, TMax, hIn);
% 
% j = 8;
% [solDDE2] = ddef4(F, phi, tau, j, TMin, TMax, hIn);
% 
% j = 14;
% [solDDE3] = ddef4(F, phi, tau, j, TMin, TMax, hIn);
% 
% Points = [left:hIn:right]';
% 
% %% Figures
% Fig1 = figure(1);
% g1 = plot(solDDE1.x,solDDE1.y(:,1),'LineWidth',1.75,'Color', [171,217,233]/255,'LineStyle','-');
% hold on 
% g2 = plot(solDDE2.x,solDDE2.y(:,1),'LineWidth',1.75,'Color', [118,42,131]/255,'LineStyle','--');
% hold on 
% g3 = plot(solDDE3.x,solDDE3.y(:,1),'LineWidth',1.75,'Color', [90,174,97]/255,'LineStyle',':');
% hold on 
% ylabel('y(t)','FontSize',15); % ,'Interpreter','latex','FontSize',15)
% xlabel('Time','FontSize',15)
% legend([g1 g2 g3],'j = 3','j = 8', 'j = 14', 'Location','NorthWest') 
