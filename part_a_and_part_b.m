 clear; close all;

% General parameters
a = 9;
b = 9;

% Part A & B parameters
time_step_no = 10;

% Part C parameters
time_step_no2 = 100000;
end_time2 = 60;

% Part D parameters
f = 8.65;
initial_u = a/b;
initial_v = a/(100*b);
time_step_no3 = 100000;
end_time3 = 200;


% Part A
% Analytical form
c = log(a/(10*b))-log(9*a^2/10);
u = @(t) ((a^2)*exp(a*t+c))./(1+a*b*exp(a*t+c));

% Analytical plot
t = linspace(0, 5/a, time_step_no);

figure
plot(t,u(t))
xlabel("Time (t)")
ylabel("Population (u)")
title("Small Fish Population","FontWeight","bold")
xlim([0,5/a])
ylim([0,1])

% Part B
% Initializing and initial conditions
u_approx = zeros([1 time_step_no]);
u_approx(1) = a/(10*b);

tau = (5/a)/(time_step_no - 1);
timeline = linspace(0, 5/a, time_step_no);

du_dt = @(u, t) a*u - b*u.^2; 


% Forward Euler Method
for i = 1:time_step_no - 1
    u_approx(i+1) =  u_approx(i) + tau*du_dt(u_approx(i), timeline(i));
end


% Modified Euler Method

predictor = zeros([1 time_step_no]);
u_high = zeros([1 time_step_no]);
u_high(1) = a/(10*b);


for i = 1:time_step_no - 1
    % Forward Euler predictor step
    predictor =  u_high(i) + tau*du_dt(u_high(i), timeline(i));

    % Trapezoidal rule with predictor
    u_high(i+1) = u_high(i) + (tau/2)*(du_dt(u_high(i), timeline(i)) + du_dt(predictor, timeline(i+1)));
end

% Plotting both approximations
hold on;
plot(timeline, u_approx, Marker="+")
plot(timeline, u_high, Marker="+")
legend("Actual", "Forward Euler", "Higher Order", "Location","best")

% Forward Euler error
error_forward = abs(u(t) - u_approx);
error_forward = error_forward'

% Higher Order error
error_higher = abs(u(t) - u_high);
error_higher = error_higher'

% Part C
% System of ODEs
du_dt2 = @(u, v, t) a*u - b*u^2 - (a + b)*u*v;
dv_dt2 = @(u, v, t) (a + b)*(u*v - v/30);

% Initialise
u2 = zeros([1 time_step_no2]);
v2 = zeros([1 time_step_no2]);

u2(1) = a/b;
v2(1) = a/(100*b);

tau2 = end_time2/(time_step_no2 - 1);
timeline2 = linspace(0, end_time2, time_step_no2);

% Modified Euler Method
for i = 1:time_step_no2 - 1
    u_predictor2 = u2(i) + tau2*du_dt2(u2(i), v2(i), timeline2(i));
    v_predictor2 = v2(i) + tau2*dv_dt2(u2(i), v2(i), timeline2(i));

    u2(i+1) = u2(i) + (tau2/2)*(du_dt2(u2(i), v2(i), timeline2(i)) + du_dt2(u_predictor2, v_predictor2, timeline2(i + 1)));
    v2(i+1) = v2(i) + (tau2/2)*(dv_dt2(u2(i), v2(i), timeline2(i)) + dv_dt2(u_predictor2, v_predictor2, timeline2(i + 1)));
end

% Plotting big and small fish population
figure;
hold on
plot(timeline2, u2)
plot(timeline2, v2)
legend("Small Fish", "Big Fish")
xlabel("Time (t)")
ylabel("Population (u)")
title("Lake Population", "FontWeight","bold")

% Part D
du_dt3 = @(u, v, t) a*u - b*u^2 - (a + b)*u*v - f*u; 
dv_dt3 = @(u, v, t) (a + b)*(u*v - v/30);

u3 = zeros([1 time_step_no3]);
v3 = zeros([1 time_step_no3]);

u3(1) = initial_u;
v3(1) = initial_v;

tau3 = end_time3/(time_step_no3 - 1);
timeline3 = linspace(0, end_time3, time_step_no3);

for i = 1:time_step_no2 - 1
    u_predictor3 = u3(i) + tau3*du_dt3(u3(i), v3(i), timeline3(i));
    v_predictor3 = v3(i) + tau3*dv_dt3(u3(i), v3(i), timeline3(i));

    u3(i+1) = u3(i) + (tau3/2)*(du_dt3(u3(i), v3(i), timeline3(i)) + du_dt3(u_predictor3, v_predictor3, timeline3(i + 1)));
    v3(i+1) = v3(i) + (tau3/2)*(dv_dt3(u3(i), v3(i), timeline3(i)) + dv_dt3(u_predictor3, v_predictor3, timeline3(i + 1)));
end

figure;
hold on
plot(timeline3, u3)
plot(timeline3, v3)
legend("Small Fish", "Big Fish")
xlabel("Time (t)")
ylabel("Population (u)")
title("Lake Population with Fishing", "FontWeight","bold")
ylim([0 0.1])


% smallfishnumber = 1;
% while abs(smallfishnumber) > 0.0001
%     
%     j = (g+h)/2;
%     if solve(g)*solve(j) < 0
%         h = j;     
%     elseif solve(h)*solve(j) < 0
%         g = j;
%     end
% end
% function u = myfunc(a,b,t)
%     c = log(a/(10*b))-log(9*a^2/10);
%     u = ((a^2)*exp(a*t+c))/(1+a*b*exp(a*t+c));
% end