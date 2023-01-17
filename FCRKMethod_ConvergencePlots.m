%% Gamma distributed DDE simulation
% Thi simulation produces the convergence plots for the FCRK method and reference solutions on the log-log scale. For
% the linear and non-linear tests, we must have j as an integer to produce
% a reference solution through the linear chain technique.

close all
clear all
% Time domain and stepsize for the DDE solver
TMin = 0;
TMax = 10;
hIn = 1e-2;

%% For the linearised differential equation
%history function
% phi = @(t) 1;

% Parameters for the distribution
% tau =  1; % 1.25;
% j = 7.0; % 7.5; % 3.85;
 
% % The RHS of the DE is in odefun where d is the delayed argument via the convolution integral
% F = @(x,d) 0.8*x-1.1*d;

%% Nonlinear comparison
tau =  2.25; % 1.25;
j = 14.0; % 3.85; 

PA.K = 2; % Also in compareLInfinity script
F = @(x,d) x-x*d/PA.K;  
phi = @(t) 0.5; 

% [sol, mySol, times, ERRGlob] = compareLinfty(F, phi, tau, j, TMin, TMax, 1e-1, 6, 4, 1, 0, 2) 
%% Comparison for linearization

tau =  4.25; % 1.25;
j = 2.25; % 7.5; % 3.85;

a = j/tau;
beta = 0.5; 
Lambda = (beta*a^j)^(1/(j+1))-a; 
phi = @(t) 1*exp(Lambda*t);

% % The RHS of the DE is in odefun  where d is the delayed argument via the convolution integral
F = @(x,d) -(j/tau).*x + beta*d;

exactSol = @(t) 1*exp(Lambda*t);

[sol, mySol, times, ERRGlob] = compareLinfty(F, phi, tau, j, TMin, TMax, 1e-1, 6, 4, 1, 1, exactSol)


%% Exact solution of linear differential equation
% tau =  1;  
% j = 1.0;  
% 
% % The RHS of the DE is in odefun  where d is the delayed argument via the convolution integral
% F = @(x,d) 0.8*x-1.1*d;
% 
% exactSol = @(t) exp(-t./10).*(cos(sqrt(29).*t./10)- (2/sqrt(29)).*sin(sqrt(29).*t./10) );
% 
%  [sol, mySol, times, ERRGlob] = compareLinfty(F, phi, tau, j, TMin, TMax, 1e-1, 6, 4, 1, 1, exactSol)

function [sol, mySol, times, ERRGlob] = compareLinfty(F, phi, tau, j, TMin, TMax, maxH, maxSteps, order, probNumber, exact, exactsol)
% This is a function that generates a convergence error plot for FCRK  solvers with fixed stepsize. 
% It compares the FCRK method against an exact one, or one found numerically with a stock
% MATLAB ode solver.
% It then saves a .pdf file of the error plot with a file name
% containig the test problem number the user decides.

% "exact" should be either '0' or '1', where '1' indicates we have the exact solution for the IVP, and '0' means we use the ODE solver.

% In the case 'exact == 0', we use ode45 to approximate the solution to the erlang IVP, and use this as the solution we compare our solvers' results to. The field 
% 'exactSol' may be filled with anything in this case--I like to use an integer !=1, but I'm pretty sure anything will work, because that value is never referenced unless 'exactSol==1'.

% In the case 'exact ==1', we feed in the exact solution as an anonymous function (i.e. exactSol = @(t) ...)
% Note that the user has to manually change the TMax hand side funcntion 'f' in the file "odefun.m" that ode45 uses, for each different IVP.

% There are error checks that could be done but aren't, like 'number of steps < 0' or 'order < 0'.

    hist = strrep(char(phi), '@(t)', '');
    rhs = strrep(char(F), '@(x,d)', '');
    %phi is a function handle, multiply it against the gamma distribution function
    if exact == 0
        a = j/tau;
        y_0 = zeros(j+1, 1);
        y_0(1) = phi(0);
        t = sym('t');
        mySol=0;
        PA.K = 2;
        for i = 2 : j+1
            %the contents of this for loop are for finding the vector of initial
            %values for ode23 or ode45, as the case may be
            %        g = @(t) exp(-a.*t) .* t.^(i-2);
            %        h = @(t) g(t).*phi(-t);
            %        initialVal = integral(h, 0, Inf,'RelTol',1e-15);
            g(t) = exp(-a*t) .* t^(i-2) ;
            h(t) = g.*phi(-t);
            initialVal = eval(int(h, t, 0, Inf));
            y_0(i) = (a.^(i-1))./gamma(i-1) .*initialVal/a;
        end
    end
    %these next two lines compute the ode23 and RK2Solver solutions (if we
    %don't provide an exact solution. we don't bother using ode45 if we are
    %only testing a 2nd order or lower fixed stepsize method
    options = odeset('RelTol',1e-15, 'AbsTol', 1e-13, 'MaxStep', 0.001);
    if order <= 2 && exact == 0
        %=sol = ode23(@(t,y) odefun(t, y, j, tau), [TMin, TMax], y_0, options);
        K= PA.K;
        sol = ode23(@(t,y) odefunNL(t, y, j, tau,K), [TMin, TMax], y_0, options);
    elseif order >2  && exact == 0
        % sol = ode45(@(t,y) odefun(t, y, j, tau), [TMin, TMax], y_0, options);
        K= PA.K;
        sol = ode45(@(t,y) odefunNL(t, y, j, tau,K), [TMin, TMax], y_0, options);
    end
    %----------------------------------------------------------------------
    %This section is strictly for visualizing the graph of the solution
    %itself compared to the actual solution 
%     Y = RK2Solver_v2(F, phi, tau, j, TMin, TMax, steps);
%     
%     queryPoints = linspace(TMin, TMax, numQueryPoints);
%     
%     ode23Eval = deval(sol, queryPoints, 1)';
%     RKEval = evalSol(Y, queryPoints');
%     
%     ERR = zeros(10, 2); %premake the vector which will contain the error with step size
%---------------------------------------------------------------------------     '
    times = zeros(maxSteps + 1, 2); % can comment this out, I had this here to see how script was taking on each loop iteration.
    steps = ceil((TMax - TMin)/maxH);
   % steps = 2*163840;
    for i = 0 : maxSteps
        if i ~=0
            steps = 2*steps;
        end
        disp(['steps: ', num2str(steps)]) %show number of steps for current loop iterate at runtime (I like to know it's working)
        h = (TMax - TMin)/(steps);
        tic;
%--------------------------------------------------------------------------
        % get my numerical solution with the number of steps specified by
        % the loop iteration, for the convergence order the caller specifies
        switch order
            case 1
                mySol = ddef1(F, phi, tau, j, TMin, TMax, h);
            case 2
                mySol = ddef2(F, phi, tau, j, TMin, TMax, h);
            case 3
                mySol = ddef3(F, phi, tau, j, TMin, TMax, h);
            case 4
                mySol = ddef4(F, phi, tau, j, TMin, TMax, h);
        end
        

 
        if exact == 0 %if we don't have an exact solution, compare with the numerical one we generated earlier (ode23 or ode45's solution is used)
            diff = abs(mySol.y(:, 1) - deval(sol, (mySol.x)', 1)');
        else
            sol = exactsol(mySol.x);
            diff = abs(mySol.y(:, 1) - exactsol(mySol.x));
        end
        MAX = max(diff);
        ERRGlob(i + 1, 1) = log(h)./log(10); %fill first column with the step size
        ERRGlob(i + 1, 2) = max(log(MAX)./log(10),log(eps)./log(10)); %fill the second column with the error
    end

    FIG = figure();
    LinearFitIndex = find(ERRGlob(:,2) > log(eps)./log(10)); % Only fit the line to points with L1 error greater than machine precision
    p1 = polyfit(ERRGlob(LinearFitIndex, 1),ERRGlob(LinearFitIndex,2), 1);
    
    fit1 = @(x) p1(2) + p1(1).*x;
    
    fplot(fit1, [ERRGlob(size(ERRGlob, 1), 1), ERRGlob(1, 1)] , 'b','LineWidth',1.25, 'DisplayName', ['Global Error Slope = ', num2str(p1(1))]);
    hold on; 
    scatter(ERRGlob(:, 1),ERRGlob(:,2), 25,'r','O','DisplayName', 'Global Error');%'ode', num2str(order), num2str(order+1)]);
%     scatter(ERRLoc (:, 1), ERRLoc(:, 2), '.', 'DisplayName', 'Local Error');

%     fplot(fit2, [ERRGlob(size(ERRGlob, 1), 1), ERRGlob(1, 1)], 'DisplayName', ['Local Error Slope = ', num2str(p2(1))]);
     [hleg] = legend('Location','NorthWest');
    % legend('Location','best','NumColumns',1,'FontSize',5);
    hold off;
    %scatter(queryPoints', ode23Eval(:,1), '.','DisplayName', 'ode23');
%     hold on 
%     plot(queryPoints', RKEval(:,1),'DisplayName', 'RK2Solver');
%     legend;
%     hold off 
    
    %title(['RK', num2str(order), ' ', 'Solver and ode', num2str(order), num2str(order+1), ' error for history function ', hist, ', RHS: ', rhs, ', j= ', num2str(j), ', tau= ', num2str(tau)])
   % title(['Order ', num2str(order), ' ', 'Solver Error' ]);
    xlabel('log_{10}(h)');
    ylabel('log_{10}(max |y(t) - u_{h}(t)|)');
    fileName2 = ['testProblem', num2str(probNumber), 'order', num2str(order), '.pdf'];

end
