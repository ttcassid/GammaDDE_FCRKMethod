function dydt = odefunNL(t, y, j, tau,K)
%This is an auxiliary function used by the solvers odexx supplied by
%MATLAB. Use this file to define erlang distributed IVP's. When changing
%the right-hand-side function f in defining a new IVP, write in the line
%marked with a (*) the function f defined in terms of x(t) and the delay
%term d. For example, if in the IVP f = x(t) + \int_{-\infty}^t x(s) *
%g_a^j(t-s) ds, we write "x + d".
    a = j/tau;
    x = y(1);
    d = a*y(j+1);
    dydt = [    x-x*d/K;    % (*)
            y(1) - a*y(2)
            ];
    for i = 3 : j+1
        dydt = [dydt ; a*(y(i-1) - y(i))];
    end
end