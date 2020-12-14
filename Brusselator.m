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

backEul = backEuler(xdt, ydt, 1, 1, 1, 3, 0.1, 30);
frwdEul = forwdEul(1, 1, 1, 3, 0.1, 30);

% plotting the solution of Brusselator problem using the Forward Euler 
% method 

% function which implements the Backward Euler method
function result = backEuler(fx, fy, X_init, Y_init, A, B, t_delta, num_iter)
data = zeros(num_iter, 3);
z_next = zeros(2, 1);
z_prev = zeros(2, 1);
F = zeros(2, 1);
z_prev(1, 1) = X_init;
z_prev(2, 1) = Y_init;
for i = 1:num_iter
   data(i, 1) = i; % the first column of data is the time step 
   % using newton's method to find the value of f_next (Z^{n+1})
   z_next = newton(z_prev(1, 1), z_prev(2, 1), A, B, num_iter);
   % f_next = (newton(X_init, Y_init, A, B, num_iter) - fz_n) / t_delta;
   % F(1, 1) = fx(z_newton(1, 1), z_newton(2, 1));
   % F(2, 1) = fy(z_newton(1, 1), z_newton(2, 1));
   % z_next = z_prev + t_delta * F;
   z_prev = z_next;
%    disp("z_next is");
%    disp(z_next);
   data(i, 2) = z_next(1, 1); % the second column of data is the value of x at time i 
   data(i, 3) = z_next(2, 1); % the third column of data is the value of y at time i 
end 
result = z_next;
end

function result = newton(X_init, Y_init, A, B, num_iter)
X = X_init;
Y = Y_init;
w_prev = zeros(2, 1);
w_next = zeros(2, 1);
w_prev(1, 1) = X_init;
w_prev(2, 1) = Y_init;
for i = 1:num_iter
    % w_next = w_prev - jacobi(X, Y, B)\F(X, Y, A, B);
    % w_next = w_prev - inv(jacobi(X, Y, B)) .* F(X, Y, A, B);
    w_prev = w_next;
    X = w_prev(1, 1);
    Y = w_prev(2, 1);
end
result = w_prev;
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
disp("The Brusselator function is");
disp(df);
end 


function h = testFunc(x, y)
h = zeros(2, 1);
h(1) = 2 * y - x;
h(2) = x;
end 


% jacobi function computes the Jacobian matrix given parameters 
% X, Y, B, and num_iter
function J = jacobi(X, Y, B)
J = zeros(2, 2);
    J(1, 1) = 2*X*Y - B - 1;
    J(1, 2) = X^2;
    J(2, 1) = B - 2*X*Y;
    J(2, 2) = - X^2;
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

