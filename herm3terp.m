function hermCoeffs = herm3terp(t_n, t_next, y_n, y_next, yp_n, yp_next)
%this is a helper function to generate the hermite polynomial's
%coefficients for the solution to the RK method over th interval

    tspan = t_next - t_n;
    slope = (y_next - y_n)/(tspan);

    %find the divided differences first
    a = y_n;
    b = yp_n;
    c = (slope - yp_n)./tspan;
    d = (yp_next - 2 .*slope + yp_n)/(tspan.*tspan);

    %premultiply the time points
    z01 = t_n .* t_n;
    z02 = t_n .* t_next;
    z012 = z01 .* t_next;

    %build the coefficients of the hermite interpolant
    a_0 = a - d .* z012 - b .* t_n + c.* z01;
    a_1 = b - 2.*c.*t_n + d.*z01 + 2.*d.*z02;
    a_2 = c - 2.*d.*t_n - d.*t_next;
    a_3 = d;

    hermCoeffs = [a_0,a_1,a_2,a_3];

end