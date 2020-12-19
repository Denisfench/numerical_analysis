close all;
% Part c: solving the Brusselator model with particular parameters
A = 1;
B = 3;
X_init = 1;
Y_init = 1;
t_delta = 0.1;
T = 30;

% test function to test Newton's method 
t = @(x, y) 2 * x^2 + 8 * y^2 - 8 - y + sqrt(3)/2 * x^2;

% constructing the differential equations 
xdt = @(X, Y) A + X^2 * Y - B * X - X;
ydt = @(X, Y) B * X - X^2 * Y;

% creating a test system of equations 
dxt = @(x, y) 2 * y - x;
dyt = @(x, y) x;
% forwardEulerTest = forwdEul(dxt, dyt, 1, 0, 1, 50);

backEul = backEuler(1, 1, 1, 3, 0.01, 30);
% frwdEul = forwdEul(1, 1, 1, 3, 0.1, 30);

% the backEuler function implements the Backward Euler method for Bruesselator
% problem, taking the initial concentrations X_0 and Y_0, constants A, B, 
% and the iteration time step df. ft is the final time T, T = 30.
function z = backEuler(X_init, Y_init, A, B, dt, ft)
num_iter = ft/dt; % the number of iterations for the backwardEler
% z = [time step, X value, Y value]
z = zeros(num_iter, 3);
z(1, 1) = 0; % initialize time = 0
z(1, 2) = X_init; % X initial 
z(1, 3) = Y_init; % Y initial 
% an actual Backward Euler loop, up to the final time T = 30 
for i = 1:num_iter - dt 
     % using Newton's method to solve the nonlinear system, in order to 
     % find the value of (Z^{n+1})
     newtn = newton(z(i, 2), z(i, 3), A, B, 10, dt);
     z(i+1, 2) = newtn(1, 1); % adding the value of x at the next time step 
     z(i+1, 3) = newtn(1, 2); % adding the value of y at the next time step
     z(i+1, 1) = z(i, 1) + dt; % adding the next time step 
end 
% plotting the time step t vs the backward Euler approximation at the 
% corresponding step t 
figure("Name", "Backward Euler: plot of time vs X and Y");
% plotting time vs X
plot(z(:, 1), z(:, 2), 'color', 'red', 'LineWidth', 2.0);
xlabel("T");
ylabel("X, Y");
hold on;
% plotting time vs Y
plot(z(:, 1), z(:, 3), 'color', 'blue', 'LineWidth', 2.0);
legend("species X", "species Y");
hold off;
end

% the newton function solves the nonlinear equation for the Backward Euler
% method, it takes X_init and Y_init as parameters, 
% which are the apporximations to the values of X and Y at a given time step 
% at which newton is called, A, B are constants for our ODE, t is the time 
% increment
function sol = newton(X_init, Y_init, A, B, num_iter, t)
% w = [X_{n+1}^{k}, X_{n+1}^{k}]
w = zeros(num_iter, 2);
% entry 1 corresponds to x value, entry 2 corresponds to y value 
w(1, 1) = X_init; % initialize the X guess for the Newton's method, which 
% corresponds to the value of X retuned by the backward Euler method at the 
% previous iteration 
w(1, 2) = Y_init; % initialize the Y guess 
% loop num_iter times 
for i = 1:num_iter-1
    % evaluating the system of ODEs at X_{n+1}^{k}, Y_{n+1}^{k}
    F = bruesselator(w(i, 1), w(i, 2), A, B);
    % evaluating the backward Euler nonlinear equation at X_{n+1}^{k}, Y_{n+1}^{k}
    G(1, 1) = w(i, 1) - X_init - t * F(1, 1);
    G(2, 1) = w(i, 2) - Y_init - t * F(2, 1);
    % performing the actual Newton's solve  
    temp = transpose(w(i, :)) - inv(jacobi(w(i, 1), w(i, 2), B, t)) * G;
    w(i+1, 1) = temp(1, 1); % adding the value of X at the next iteration
    w(i+1, 2) = temp (2, 1); % adding the value of Y at the next iteration
end
sol = w(num_iter, :);
end 


% defining the function dX/dt of the first ODE
function f_xdt = f_xdt(x, y, a, b)
f_xdt = a + x^2 * y - b * x - x;
end 


% defining function dY/dt of the second ODE
function f_ydt = f_ydt(x, y, b)
f_ydt = b * x - x^2 * y;
end 


% using dX/dt and dY/dt to construct the sytem of ODEs
function F = bruesselator(x, y, a, b)
F = zeros(2, 1);
F(1) = f_xdt(x, y, a, b);
F(2) = f_ydt(x, y, b);
end 


% The function below implements the Jacobian matrix for the Newton's step
% the Jacobian is for the fucntion G(Z^{n+1})
function J = jacobi(X, Y, B, t)
    J = zeros(2, 2);
    J(1, 1) = 1 - 2 * X * Y * t + B * t + t;
    J(1, 2) = - X^2 * t;
    J(2, 1) = - B * t + 2 * X * Y * t;
    J(2, 2) = 1 + X^2 * t;
end

% Forward Euler method implementation 
function z = forwdEul(X_init, Y_init, A, B, h, num_iter)
% z is a 3D vector first column of which corresponds to the time step being
% taken, the second column corresponds to the value of x, and the third
% column corresponds to the value of y at a particular time step 
temp = zeros(2, 1);
z = zeros(num_iter, 3);
F = zeros(2, 1);
z(1, 1) = 1;
z(1, 2) = X_init;
z(1, 3) = Y_init;
F = bruesselator(z(1, 2), z(1, 3), A, B);
for i = 1:num_iter-1
     % y{i+1} = y_{i} + h * F(x{i}, y{i});
     temp = z(i, 3) + h * F; % F is 2D 
     % update the value of x 
     z(i+1, 2) = temp(1, 1);
     % update the value of y 
     z(i+1, 3) = temp(2, 1);
     % set the value of time step T
     z(i+1, 1) = i+1;
     F = bruesselator(z(i+1, 2), z(i+1, 3), A, B);
end 
% plotting the time T vs forward Euler approximation
figure("Name", "Forward Euler: plot of X and Y vs time");
% plotting x vs time 
plot(z(:, 1), z(:, 2), 'color', 'red', 'LineWidth', 2.0);
xlabel("T");
ylabel("X, Y");
hold on;
% plotting y vs time 
plot(z(:, 1), z(:, 3), 'color', 'blue', 'LineWidth', 2.0);
legend("species X", "species Y");
hold off;
end

