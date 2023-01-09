function sol = ddef4(F, phi, tau, j, left, right, hIn)
%   WHAT'S NEW IN THIS VERSION
%   changed getMesh to version 4

%   A program to solve the initial value problem y'(t) = F(y(t), D(t)) where
%   D is a distributed delay. F must be entered as an anonymous function F =
%   @(x, D). phi is the history function for the IVP, also an anonymous
%   function. left and right are the endpoints of the solution interval for
%   the IVP, and H is the step size the solver will use
%   to find the solution y(t) on [left, right].
    
%--------------------------------------------------------------------------
%    INPUT VALIDATION

if ~isa(F,'function_handle')
    ME = MException('InputException:TypeError',...
        'F is not defined as a function handle.');
    throw(ME)
end
if ~isa(phi,'function_handle')
    ME = MException('InputException:TypeError',...
        'phi is not defined as a function handle.');
    throw(ME)
end
if j <= 0
    ME = MException('InputException:OutOfBounds', ...
        'Variable j is <=0');
    throw(ME)
end
if tau <= 0
    ME = MException('InputException:OutOfBounds', ...
        'Variable tau is <=0');
    throw(ME)
end
if right <= left
    ME = MException('InputException:OutOfBounds', ...
        'Interval [right, left] not defined or degenerate. Ensure right > left.');
    throw(ME)
end
if hIn <= 0
    ME = MException('InputException:OutOfBounds', ...
        'Enter a step size h > 0');
    throw(ME)
end
if hIn > right - left
    ME = MException('InputException:OutOfBounds', ...
        'Stepsize h > right - left; stepsize larger than interval');
    throw(ME)
end
%--------------------------------------------------------------------------

%   Defining constants to be used throughout the script
    numSteps = floor((right - left)./hIn) + 1; %find the number of steps 
    a = j/tau;
    h = (right - left) /numSteps; %define the step size for IVP solution
    k = 1 + 1/j;
    beta = 1 + 5/j;
    A = (j+1)/a^(1/beta); % to be passed to the simpson's rule function
    g= @(t)((beta.*a.^j.* (k.*j./a.^(1./beta)).^(beta.*j))./gamma(j)).* exp(-a.*(-k.*j.*log(t)./a.^(1./beta)).^beta).*(-log(t)).^(beta .*j-1).*(1./t);
    
%-------------------------------------------------------------------------

%   global variables pertaining to numerical quadrature
    intCoeffs = [2, -1, 2, 2, -1, 2, 2, -1, 2]; %this will eventually be a vector of [2, -1, 2] growing as necessary
    K = [0, 0, 0, 0, 0, 0]; %this is a global variable containing the stage values K_i
    hQuadrature = 0.5 .* h; %this is the step size which will be used for numerical integration 
    t_right = [0,0];
    specialPoint= [0,0];
    reuseIndex = [0, 0];
    intMeshSize = [0,0];
    evalFirstPart1 = [];
    I_numPartial = [0, 0];
    I_Stage = [0,0,0,0,0,0];
%-----------------------------------------------------------------------
    
%   create space to store the solution
    SOL = zeros(numSteps+1, 6); %preallocate solution for speed. 
    SOL(1, 1) = left;
    SOL(1, 2) = phi(left); %set the value of history function at t_0
        
%-----------------------------------------------------------------------
    
%   find K(1) for the first time
    t_n = left;
    y_n = SOL(1, 2);
    quadMesh1 = getMesh_v4(hQuadrature, 1);
    transformedMesh1 = t_n - (-A.*log(quadMesh1.mesh)).^(beta);
    intMeshSize(1) = size(transformedMesh1, 1);
    gMesh1 = g(quadMesh1.mesh); 
    gMeshH1 = multiplyByH(gMesh1, quadMesh1.h, quadMesh1.hTemp1, quadMesh1.hTemp2, quadMesh1.index, intMeshSize(1));
    [evalFirstPart1, reuseIndex(1)] = assignFirstPart(SOL, transformedMesh1, intMeshSize(1), left, t_n, phi);
    I_numPartial(1) = integrate(1, gMeshH1, evalFirstPart1);
    K(1) = F(y_n, I_numPartial(1));
    yp_n = K(1);
    
%   main loop
    for i = 1 : numSteps
% set t_right
        t_right(1) = t_n + (1/2).*h;
        t_right(2) = t_n + h;
%       get new meshes for numerical quadrature
        specialPoint(1) = exp((-1/A).*(t_right(1) - left).^(1/beta));
        specialPoint(2) = exp((-1/A).*(t_right(2) - left).^(1/beta));       
        quadMesh1 = getMesh_v4(hQuadrature, specialPoint(1)); %generate the mesh which will be used for the integration up to t_right
        quadMesh2 = getMesh_v4(hQuadrature, specialPoint(2)); %generate the mesh which will be used for the integration up to t_right
        gMesh1 = g(quadMesh1.mesh); 
        gMesh2 = g(quadMesh2.mesh); 
%       transform the meshes to be over (-\infty, t_right)
        transformedMesh1 = t_right(1) - (-A.*log(quadMesh1.mesh)).^(beta);
        transformedMesh2 = t_right(2) - (-A.*log(quadMesh2.mesh)).^(beta);

%       Assign function values to the mesh points which fall in (-\infty, t_n)
%       first call 'assignFirstPart' which will assign all function
%       values on (-\infty, t_n)
        intMeshSize(1) = size(transformedMesh1, 1);
        intMeshSize(2) = size(transformedMesh2, 1);
        [evalFirstPart1, reuseIndex(1)] = assignFirstPart(SOL, transformedMesh1, intMeshSize(1), left, t_n, phi, h); %reuseIndex is the index of the first point which must be evaluated at time > t_n
        [evalFirstPart2, reuseIndex(2)] = assignFirstPart(SOL, transformedMesh2, intMeshSize(2), left, t_n, phi, h);
%       once reuseIndex is known, we can pre-multiply each mesh point by 
%       its corresponding hQuadrature (will only be different for the
%       interval affected by the specialPoint
        gMeshH1 = multiplyByH(gMesh1, quadMesh1.h, quadMesh1.hTemp1, quadMesh1.hTemp2, quadMesh1.index, intMeshSize(1));
        gMeshH2 = multiplyByH(gMesh2, quadMesh2.h, quadMesh2.hTemp1, quadMesh2.hTemp2, quadMesh2.index, intMeshSize(2));

%       for t_n + h/2 and t_n + h, integrate the solution up to time t_n
        I_numPartial(1) = integrate(1, gMeshH1, evalFirstPart1);
        I_numPartial(2) = integrate(1, gMeshH2, evalFirstPart2);
     
%       To evaluate the transformed mesh at the stage interopolants, need
%       to convert from t to alpha = (t - t_n)/h \in [0, c_i], for each
%       element of t_right. Efficient to compute these vectors here once per RK step instead of every time
%       assignSecondPart() is called.

        alph1 = (transformedMesh1(reuseIndex(1) : intMeshSize(1)) - t_n);
        alph1sq = alph1 .* alph1;
        alph1cub = alph1 .* alph1sq;
        alph2 = (transformedMesh2(reuseIndex(2) : intMeshSize(2)) - t_n);
        alph2sq = alph2.*alph2;
        alph2cub = alph2.*alph2sq;
        
%       then call 'assignSecondPart' for subsequent stages, integrate the
%       portion from t_n to t_n + h or t_n + h/2 as the case may be, and
%       compute K(i)

        stageEval1 = assignSecondPart (alph2, alph2sq, alph2cub, intMeshSize(2), reuseIndex(2), y_n, 1, K, h);
        I_Stage(1) = integrate(reuseIndex(2), gMeshH2, stageEval1);
        K(2) = F(y_n + h.*K(1), I_numPartial(2) + I_Stage(1));

        stageEval2 = assignSecondPart(alph1, alph1sq, alph1cub, intMeshSize(1), reuseIndex(1), y_n, 2, K, h);
        I_Stage(2) = integrate(reuseIndex(1), gMeshH1, stageEval2);
        K(3) = F(y_n + h.*(K(1).*(3./8) + K(2)./8), I_numPartial(1) + I_Stage(2));
        
        stageEval3 = assignSecondPart(alph2, alph2sq, alph2cub, intMeshSize(2), reuseIndex(2), y_n, 3, K, h);
        I_Stage(3) = integrate(reuseIndex(2), gMeshH2, stageEval3);
        K(4) = F(y_n + h.*(K(1)./2 + K(2)./2), I_numPartial(2) + I_Stage(3));
        
        stageEval4 = assignSecondPart(alph1, alph1sq, alph1cub, intMeshSize(1), reuseIndex(1), y_n, 4, K, h);
        I_Stage(4) = integrate(reuseIndex(1), gMeshH1, stageEval4);
        K(5) = F(y_n + h .*( (5./24).*K(1) + K(3)./3 - K(4)./24 ) , I_numPartial(1) + I_Stage(4));
        
        stageEval5 = assignSecondPart(alph2, alph2sq, alph2cub, intMeshSize(2), reuseIndex(2), y_n, 5, K, h);
        I_Stage(5) = integrate(reuseIndex(2), gMeshH2, stageEval5);
        K(6) = F( y_n + h.* ( K(1)./6 + (2./3).*K(3) + K(4)./6 ) , I_numPartial(2) + I_Stage(5) );
        
%       Now that the K_i are all known, we find y_next
        y_next = y_n + h .* ( K(1)./6 + (2./3).*K(5) + K(6)./6 );
        
%       evaluate K(1) of the endpoint t_next, but now with the final stage interpolant given by
%       the RK method as the function definition
        stageEval6 = assignSecondPart(alph2, alph2sq, alph2cub, intMeshSize(2), reuseIndex(2), y_n, 5, K, h);
        I_Stage(6) = integrate(reuseIndex(2), gMeshH2, stageEval6);
        K(1) = F(y_next, I_numPartial(2) + I_Stage(6));
        yp_next = K(1);
        
%       Update the matrix of solution information
        SOL(i+1, 1) = t_right(2);
        SOL(i+1, 2) = y_next;
%       Load the coefficients of the interpolant defined on [t_n, t_next] into the solution matrix        
        SOL(i, 3:6) = herm3terp(t_n, t_right(2), y_n, y_next, yp_n, yp_next); 
%       advance t and y
        t_n = t_right(2);
        y_n = y_next;
        yp_n = yp_next;
    end

    function I_num = integrate(leastIndex, gMeshH, assignmentValues)
%   integrates on the mesh that is supplied to it from the least index up to the highest
%   index. 
        maxIndex = leastIndex + size(assignmentValues, 1)-1;
        while size(intCoeffs, 2) < intMeshSize(1) || size(intCoeffs, 2) < intMeshSize(2)
                intCoeffs = [intCoeffs, 2, -1, 2]; %grow the intCoeffs row vector as necessary
        end
%       Perform a dot product over the mesh points to approximate the integral  needed.
%       For example, to find the integral from -\infty to t_0, only
%       dot over the first reuseIndex - 1 coefficients,
%       and the corresponding function evaluations found in
%       'assignmentValues'
        I_num = intCoeffs(1, leastIndex: maxIndex) * (assignmentValues.*gMeshH(leastIndex: maxIndex));
    end
    SOL(numSteps + 1, 1) = right; %do this or 'deval' will complain about exceeding the solution interval because of rounding errors
    sol.x = SOL(:, 1);
    sol.y = SOL(:, 2:6);
    sol.hist = phi;
    sol.interpolantOrder = 4;
end

function [evalFirstPart, reuseIndex] = assignFirstPart(SOL, transformedMesh, intMeshSize, t_0, t_n, phi, h)
%   take as input a mesh from (-\infty, t_right) and assign the correct evaluation to each point of
%   the transformed mesh falling between -\infty and t_n. It also outputs the
%   index in this mesh of the first point which is > t_n and < t_right

%   fnDefs has format [t_n, yn, ypn, a_0, a_1] where a_i is a
%   coefficient of the linear interpolant on the interval starting with t_n
%   t_left is the left endpoint of current interval we are solving in.
       
       evalFirstPart = zeros(intMeshSize, 1);
       reuseIndex = 0;
       for loopIndex = 1 : intMeshSize
            currPoint = transformedMesh(loopIndex);
            if currPoint <= t_0 %if the integration mesh point is in the history, evaluate it as such
                evalFirstPart(loopIndex) = phi(currPoint);
            elseif currPoint > t_0 && currPoint <= t_n 
%               i.e. if the transformed point falls somewhere in the domain of
%               the already computed solution
                interval = floor((currPoint - t_0)/h)+1;
                evalFirstPart(loopIndex) =  SOL(interval, 3 : 6 ) * [1; currPoint; currPoint.*currPoint; currPoint.*currPoint.*currPoint];
            elseif currPoint > t_n %apply the stage interpolant for the current interval
                   reuseIndex = loopIndex; 
                break
            end
       end
       if reuseIndex ~= 0 %only return the relevant part of EvalFirstPart
            evalFirstPart = evalFirstPart(1: reuseIndex - 1);
       end
end

function stageEval = assignSecondPart(t, t2, t3, intMeshSize, reuseIndex, y_n, stageNumber, K, h)
%   This function takes as input a mesh from (-\infty, t_right) and the index
%   of the first point which is > t_n and < t_right. Depending on the stage,
%   use a different function to evaluate these points.
    stageEval = zeros(intMeshSize - reuseIndex +1, 1); %allocate space for the stageEval vector
    switch stageNumber
        case 1
            stageEval = y_n + (K(1)).*(t);
        case {2 , 3}
            a2 = (1/2).* t2./h;
            a1 = t - a2; 
            stageEval= y_n + ( a1.*K(1)+ a2.* K(2));
        case {4 , 5}
            a1 = t - (3./(2.*h)).*t2 + (2./(3.*h^2)).*t3;
            a3 = (2./h).*t2 - (4./(3.*h^2)).*t3;
            a4 = -t2./(2.*h) + (2./(3.*h^2)).*t3;
            stageEval= y_n +(  K(1).*a1 + K(3).*a3 + K(4).*a4  );
        case 6
            a1 = t - (3./(2.*h)).*t2 + (2./(3.*h^2)).*t3;
            a5 =(2./h).*t2 - (4./(3.*h^2)).*t3;
            a6 = -t2./(2.*h) + (2./(3.*h^2)).*t3;
            stageEval = y_n + ( K(1).*a1 + K(5).*a5 + K(6).*a6 );
    end 
end

function gMeshH = multiplyByH(gMesh, hDefault, hTemp1, hTemp2, spIndex, gMeshSize)
%   This function multiplies each entry in the gMesh by its appropriate
%   step size. This is slightly complicated since the presence of the
%   special point creates six mesh points with a different spacing
    gMeshH = zeros(gMeshSize, 1);
    if spIndex < 0 % there is no special point, it's already in the mesh, and so every mesh point gets the defeault h
        gMeshH = hDefault/3 .* gMesh;
    else
        gMeshH(1 : spIndex - 1) = hDefault./3.*gMesh(1 : spIndex-1);
        gMeshH(spIndex : spIndex + 2) = hTemp1./3 .*gMesh(spIndex : spIndex + 2);
        gMeshH(spIndex+3 : spIndex + 5) = hTemp2./3 .*gMesh(spIndex+3 : spIndex + 5);
        gMeshH(spIndex + 6 : gMeshSize) = hDefault./3.*gMesh(spIndex + 6 : gMeshSize);
    end
end

function newMesh = getMesh_v4(hQuadrature, specialPoint)
    numPoints = floor((1/hQuadrature))+1; %determines a uniform mesh size which is slightly less than the one pre-calculated in the preamble of the script
    h = 1 ./ numPoints;
    %divide h by 4 once and save these values to reduce rounding error
    h_1 = h./4;
    h_2 = h./2;
    h_3 = 3*h_1;
    
    test = numPoints * specialPoint;
    spInterval = 0; %this will contain the interval of the mesh which contains the special point
    index = 0;
    if  test == floor(test) %check to see if the special point is already part of the mesh, If so, make the mesh without the special point
        mesh = zeros(3 * (numPoints), 1);
        index = -1;
        spInterval = -1;
        hTemp1 = -1;
        hTemp2 = -1;
        for i = 0 : 3: 3*(numPoints)-1
            mesh(i+1, 1) = h*(floor(i/3)) + h_1;
            mesh(i+2, 1) = h*(floor(i/3)) + h_2;
            mesh(i+3, 1) = h*(floor(i/3)) + h_3;
        end
    else 
        mesh = zeros(3*numPoints +3, 1); %change number of columns to 2 if you want to debug spacing of mesh points 
        spInterval = floor(test); %this places the special point in its proper interval
        counter = 0;
        interval = 0;
        while interval < spInterval && interval < numPoints-1
            mesh(counter+1, 1) = h * interval + h_1;
            mesh(counter+2, 1) = h * interval + h_2;
            mesh(counter+3, 1) = h * interval + h_3;
            counter = counter + 3;
            interval = interval + 1;
        end
        index = counter+1; %save the index of first meshpoint in the left 'specialPoint' interval. 
        hTemp1 = specialPoint - h*interval;
        mesh(counter+1, 1) = h * interval + hTemp1 /4;
        mesh(counter+2, 1) = h * interval + hTemp1 /2;
        mesh(counter + 3, 1) = h * interval + 3*hTemp1 / 4;
         counter = counter + 3;
        %spInterval = interval; % save the mesh interval later to add to the mesh struct in the output
        interval = interval + 1;
        
        hTemp2 = h*interval - specialPoint;
        mesh(counter+1, 1) = specialPoint + hTemp2 /4;
        mesh(counter+2, 1) = specialPoint + hTemp2 /2;
        mesh(counter + 3, 1) = specialPoint + 3*hTemp2 / 4;
        counter = counter + 3;
        
        while interval < numPoints
             mesh(counter+1, 1) = h * interval + h_1;
            mesh(counter+2, 1) = h * interval + h_2;
            mesh(counter+3, 1) = h * interval + h_3;
            counter = counter + 3;
            interval = interval + 1;
        end
    end
    newMesh.index = index;
    newMesh.numPoints = numPoints;
    newMesh.mesh = mesh;
    newMesh.interval = spInterval;
    newMesh.h = h;
    newMesh.hTemp1 = hTemp1;
    newMesh.hTemp2 = hTemp2;
end