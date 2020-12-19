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

backEul = backEuler(xdt, ydt, 1, 1, 1, 3, 0.01, 30);
% frwdEul = forwdEul(1, 1, 1, 3, 0.1, 30);

% plotting the solution of Brusselator problem using the Forward Euler 
% method 

% function which implements the Backward Euler method
function z = backEuler(fx, fy, X_init, Y_init, A, B, dt, ft)
num_iter = ft/dt; % the number of iterations is the 
z = zeros(num_iter, 3);
z(1, 1) = 0; % time why is the first time step 1, not 0? it should be 0 and
% increase in increments of dt 
z(1, 2) = X_init; % x initial 
z(1, 3) = Y_init; % y initial 
for i = 1:num_iter - dt
     % using newton's method to find the value of (z^{n+1})
     newtn = newton(z(i, 2), z(i, 3), A, B, 10, dt);
     disp(newtn);
     % temp = z(i, [2, 3]) + t_delta * bruesselator(z(i, 2), z(i, 3), A, B);
     % disp(z(i,:));
     % update the value of x 
     z(i+1, 2) = newtn(1, 1);
     % update the value of y 
     z(i+1, 3) = newtn(1, 2);
     % set the value of time step T
     z(i+1, 1) = z(i, 1) + dt;
end 
% plotting the time T vs forward Euler approximation
figure("Name", "Backward Euler: plot of X and Y vs time");
% plotting time vs x 
plot(z(:, 1), z(:, 2), 'color', 'red', 'LineWidth', 2.0);
xlabel("T");
ylabel("X, Y");
hold on;
% plotting time vs y 
plot(z(:, 1), z(:, 3), 'color', 'blue', 'LineWidth', 2.0);
legend("species X", "species Y");
hold off;
end

% testing the newton's method 
% doing 5 iterations 
function sol = newton(X_init, Y_init, A, B, num_iter, t)
temp = zeros(2, 1);
G = zeros(2, 1);
F = zeros(2, 1);
w = zeros(num_iter, 2);
% entry 1 corresponds to x value, entry 2 corresponds to y value 
w(1, 1) = X_init;
w(1, 2) = Y_init;
for i = 1:num_iter-1
    % precomputing G(Z_{n+1}^k), Z^n = [X_intit, Y_init]
    % computing F(Z_{n+1}^k)
    F = bruesselator(w(i, 1), w(i, 2), A, B);
    G(1, 1) = w(i, 1) - X_init - t * F(1, 1);
    G(2, 1) = w(i, 2) - Y_init - t * F(2, 1);
    % temp = transpose(w(i, :)) - inv(jacobi(w(i, 1), w(i, 2), B, t)) * bruesselator(w(i, 1), w(i, 2), A, B);
    temp = transpose(w(i, :)) - inv(jacobi(w(i, 1), w(i, 2), B, t)) * G;
    % disp(temp)
    w(i+1, 1) = temp(1, 1);
    w(i+1, 2) = temp (2, 1);
end
sol = w(num_iter, :);
end 

% defining function f_xdt as a function handle
function f_xdt = f_xdt(x, y, a, b)
f_xdt = a + x^2 * y - b * x - x;
end 

% defining function f_xdt as a function handle
function f_ydt = f_ydt(x, y, b)
f_ydt = b * x - x^2 * y;
end 

% Bruesselator problem definition
function df = bruesselator(x, y, a, b)
F = zeros(2, 1);
F(1) = f_xdt(x, y, a, b);
F(2) = f_ydt(x, y, b);
df = F;
%disp("The Brusselator function is");
%disp(df);
end 


function h = testFunc(x, y)
h = zeros(2, 1);
h(1) = 2 * y - x;
h(2) = x;
end 


% jacobi function computes the Jacobian matrix given parameters 
% X, Y, B, and num_iter
% I believe the Jacobian below id incorrect 
% function J = jacobi(X, Y, B)
% J = zeros(2, 2);
%     J(1, 1) = 2*X*Y - B - 1;
%     J(1, 2) = X^2;
%     J(2, 1) = B - 2*X*Y;
%     J(2, 2) = - X^2;
% end

% The function below implements the Jacobian matrix for the function 
% G(Z^{n+1})
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

