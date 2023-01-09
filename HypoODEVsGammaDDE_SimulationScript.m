%% Gamma distributed DDE simulation
close all
clear all

% Coefficients for the linearised differential equation
PA.Alpha =  0.8;  
PA.Beta =  -11/10;  

%Distribution parameters
tau =  1.0;  
j = 2.57;  

% Unstable vs stable spiral ODE not in DDE  
% PA.Alpha =   0.825; %   
% PA.Beta =   -1.175; % 
% tau =  1.0; % 1.5;
% j = 4.495;  

% Stable vs unstable spiral ODE not in DDE 
% PA.Alpha =   0.89; %   
% PA.Beta =   -1.15; %    
% tau =  1.0; % 1.5;
% j = 2.5; %   

F = @(x,d) PA.Alpha*x+PA.Beta*d;  % RHS of the linearised distributed DDE where d is the delayed argument via the convolution integral
phi = @(t) 0.1*exp(1e-1*t); 

%% Nonlinear DDE

% PA.K = 2;
% F = @(x,d) x-x*d/PA.K;  % RHS of the linearised distributed DDE where d is the delayed argument via the convolution integral
% phi = @(t) 0.5;  % history function for the IVP

% Distribution parameters 
% tau =  2.25;  
% j = 6.45;  

%% Set up the DDE solver
%Time interval (TMin & TMax endpoints+ timestep.) 
TMin = 0;
TMax = 250; % 20;
hIn = 1e-2;

% Parameters for the distribution
[solDDE] = ddef4(F, phi, tau, j, TMin, TMax, hIn);

Points = [TMin:hIn:TMax]';

%% Hypoexponential ODE 2 free rates
lambda = j/tau;
PA.N = max(ceil(j),2);

PA.TransitRate = PA.N*lambda/j; 

c2 = 1/(PA.TransitRate*lambda) ;
c1 = (PA.N-2)*( (1/PA.TransitRate)^2 - c2 );
t = sym('t');

if j == round(j) % If the underlying DDE is Erlang
    PA.XRootV1 =  PA.TransitRate;  
    PA.YRootV1 =  PA.TransitRate;  
    TransitCompartIC = zeros(1,PA.N);
    
    for ii = 1 : PA.N  %Calculate the initial conditions for transit compartments of the Erlang distribution
        g(t) = exp(-PA.TransitRate*t) .* t^(ii-1) ;
        h(t) = g.*phi(-t);
        initialVal = eval(int(h, t, 0, Inf));
        TransitCompartIC(ii) = ((PA.TransitRate.^(ii))./gamma(ii)).*initialVal/PA.TransitRate;
    end
    ICVec  = [TransitCompartIC(1:end)];
else % If the underlying DDE is Not Erlang
    G = @(x) 2 - 4*(tau/PA.N)*x + ( 4*(tau^2/PA.N^2) - 2*c2 + c1 )*x^2 ;
    PA.XRootV1 =  fzero(G,0.75*PA.TransitRate);
    PA.YRootV1 =  1./( 2*tau/PA.N - 1/PA.XRootV1 ); 
    TransitCompartIC = zeros(1,PA.N-2); 
    
    for ii = 1 : PA.N-2 %Calculate the initial conditions for transit compartments of the Hypoexponential distribution
        g(t) = exp(-PA.TransitRate*t) .* t^(ii-1) ;
        h(t) = g.*phi(-t);
        initialVal = eval(int(h, t, 0, Inf));
        TransitCompartIC(ii) = ((PA.TransitRate.^(ii))./gamma(ii)).*initialVal/PA.TransitRate;
    end
    
    gx(t) = exp(-PA.XRootV1*t).*PA.XRootV1;
    hx(t) = gx.*phi(-t)./PA.XRootV1;
    XIC = eval(int(hx, t, 0, Inf));
    
    gy(t) = exp(-PA.YRootV1*t).*PA.YRootV1;
    hy(t) = gy.*phi(-t)./PA.YRootV1;
    YIC = eval(int(hy, t, 0, Inf));
    ICVec = [TransitCompartIC(1:end),XIC,YIC];
end

if PA.YRootV1 < 0 %Warning
    disp('Warning: Transit rate not defined')
else
% Test the two moment matching    
TestFirstCond = (PA.N-2)*(1/PA.TransitRate) + 1/PA.XRootV1 + 1/PA.YRootV1 - tau;
TestSecondCond = (PA.N-2)*(1/PA.TransitRate)^2 + 1/PA.XRootV1^2 + 1/PA.YRootV1^2  - tau/lambda;

ODEApproxICV1 =  [phi(0), ICVec] ; 
totaltime = [TMin TMax];
[solODEApproxV1] = ODEApproxSolverTwoParam(totaltime,ODEApproxICV1,PA);

%% Round Erlang ODE Comparison
PA.N  =  round(j)  ;
PA.TransitRate = PA.N/tau;

TransitCompartIC = zeros(1,PA.N);

for ii = 1 : PA.N  %Calculate the initial conditions for transit compartments of the Erlang distribution
    g(t) = exp(-PA.TransitRate*t) .* t^(ii-1) ;
    h(t) = g.*phi(-t);
    initialVal = eval(int(h, t, 0, Inf));
    TransitCompartIC(ii) = ((PA.TransitRate.^(ii))./gamma(ii)).*initialVal/PA.TransitRate;
end
ICVec  = [TransitCompartIC(1:end)];

ErlangApproxICV1 =  [phi(0),ICVec]; % (phi(0)./PA.TransitRate).*ones(1,PA.N)] ;
totaltime = [TMin TMax];
[solODERoundErlang] = ODEErlangApproxSolver(totaltime,ErlangApproxICV1,PA);

%% Figures
Fig3 = figure(3);
g1 = plot(solDDE.x,solDDE.y(:,1),'LineWidth',1.75,'Color', [171,217,233]/255,'LineStyle','-');
hold on 
g2 = plot(solODEApproxV1.x,solODEApproxV1.y(1,:),'LineWidth',1.75,'Color',[239,138,98]/256,'LineStyle','--'); 
hold on  
g3 = plot(solODERoundErlang.x,solODERoundErlang.y(1,:),'LineWidth',0.75,'Color',[118,42,131]/256,'LineStyle','-'); 
hold on 
ylabel('y(t)','FontSize',15); % ,'Interpreter','latex','FontSize',15)
xlabel('Time','FontSize',15)
legend([g1 g2 g3],'DDE solution','Two Moment','Erlang','Location','SouthWest') 

Fig4 = figure(4);
RelativeError = max(abs(solDDE.y(:,1)),1e-5); 
TestSol = deval(solODEApproxV1,solDDE.x,1);
g3 = plot(solDDE.x,log( ( abs(solDDE.y(:,1)-TestSol') ) )./log(10),'LineWidth',2.25,'Color', [239,138,98]/256,'LineStyle','-');
hold on
TestSol2 = deval(solODERoundErlang,solDDE.x,1);
g3 = plot(solDDE.x,log( ( abs(solDDE.y(:,1)-TestSol2') ) )./log(10),'LineWidth',2.25,'Color', [118,42,131]/256,'LineStyle','-');
ylabel('log_{10} (|y(t)-u(t)|)','FontSize',15); % ,'Interpreter','latex','FontSize',15)
xlabel('Time','FontSize',15)
 
end



function [sol] = ODEApproxSolver(totaltime,IC,PA) %DDE model without therapy
opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep',1e-2);
sol = ode45(@ODEApproxDE,totaltime,IC,opts);
        function dydt = ODEApproxDE(t,y);
            dydt(1) =  PA.Alpha*y(1) + PA.Beta*PA.ZRoot*y(PA.N+1); % y(1) - y(1)*PA.ZRoot*y(PA.N+1)/PA.K; %   % %  + cos(y(PA.N+1)) ; % 
            if PA.N == 4
            dydt(2) = y(1) - PA.TransitRate*y(2);
            dydt(3) = PA.TransitRate*y(2) - PA.XRoot*y(3);
            dydt(4) = PA.XRoot*y(3) - PA.YRoot*y(4);
            dydt(5) = PA.YRoot*y(4) - PA.ZRoot*y(5);
            else
            dydt(2) =  y(1) - PA.TransitRate*y(2);   
            for ii = 1:PA.N-3 %y(2) is the first compartment
                dydt(2+ii) = PA.TransitRate.*(y(2+ii-1) -y(2+ii) );
            end
            dydt(PA.N-1) =  PA.TransitRate.*y(PA.N-2) - PA.XRoot.*y(PA.N-1);
            dydt(PA.N) =   PA.XRoot.*y(PA.N-1) - PA.YRoot*y(PA.N);
            dydt(PA.N+1) = PA.YRoot*y(PA.N) - PA.ZRoot*y(PA.N+1); 
            end
            dydt = dydt';
        end
end
 
function [sol] = ODEApproxSolverTwoParam(totaltime,IC,PA) %DDE model without therapy
opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep',1e-2);
sol = ode45(@ODEApproxDETwoParam,totaltime,IC,opts);
        function dydt = ODEApproxDETwoParam(t,y);
            dydt(1) = PA.Alpha*y(1) + PA.Beta*PA.YRootV1*y(PA.N+1); %  y(1) - y(1)*PA.YRootV1*y(PA.N+1)/PA.K; %    
            % 
            if PA.N == 2
            dydt(2) = y(1) - PA.XRootV1*y(2);
            dydt(3) = PA.XRootV1*y(2) - PA.YRootV1*y(3); 
            else
            dydt(2) =  y(1) - PA.TransitRate*y(2);   
            for ii = 1:PA.N-3 %y(2) is the first compartment
                dydt(2+ii) = PA.TransitRate.*(y(2+ii-1) -y(2+ii) );
            end
            dydt(PA.N) =  PA.TransitRate.*y(PA.N-1) - PA.XRootV1.*y(PA.N);
            dydt(PA.N+1) =   PA.XRootV1.*y(PA.N) - PA.YRootV1*y(PA.N+1);
            end
            dydt = dydt';
        end
end

 function [sol] = ODEErlangApproxSolver(totaltime,IC,PA) %DDE model without therapy
opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep',1e-2);
sol = ode45(@ODEApproxDETwoParam,totaltime,IC,opts);
        function dydt = ODEApproxDETwoParam(t,y);
            dydt(1) =   PA.Alpha*y(1) + PA.Beta*PA.TransitRate*y(PA.N+1); %y(1) - y(1)*PA.TransitRate.*y(PA.N+1)/PA.K;  % 
            %
            if PA.N == 1
            dydt(2) = y(1) - PA.TransitRate*y(2);  
            elseif PA.N == 2
            dydt(2) = y(1) - PA.TransitRate*y(2);
            dydt(3) =PA.TransitRate*(y(2) - y(3)); 
            else
            dydt(2) =  y(1) - PA.TransitRate*y(2);   
            for ii = 1:PA.N-1 %y(2) is the first compartment
                dydt(2+ii) = PA.TransitRate.*(y(2+ii-1) -y(2+ii) );
            end
            end
            dydt = dydt';
        end
 end