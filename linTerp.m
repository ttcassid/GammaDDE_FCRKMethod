function coeffs = linTerp(t_n, t_next, y_n, y_next)
%this is a helper function to generate the coefficients for a piecewise
%linear interpolant for the function values of the solution.

    slope = (y_next - y_n)./(t_next - t_n);
    a_0 = y_n - slope.*t_n;
    a_1 = slope;
    
    coeffs = [a_0, a_1];
end