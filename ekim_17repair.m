clc; clear all; close all;
% Test for Linearity

% Define the discrete time index n
n = -10:10;

% Define the unit step function u[n]
u = (n >= 0); % unit step function (u[n] = 1 for n >= 0), 0 otherwise

% Define the input signals in the form of a^n * u[n]
a1 = 0.9; % First base for exponential input
a2 = 0.3; % Second base for exponential input

x1 = a1.^n .* u; % First input signal: a1^n * u[n]
x2 = a2.^n .* u; % Second input signal: a2^n * u[n]

% Define the systems
for s = 1:1:4
    if s == 1
        % System 1: y[n] = x[n] + x[n+5] (linear)
        T = @(x,n) x + [x(6:end) zeros(1,5)];
    elseif s == 2
        % System 2: y[n] = a^(n-1) * x[n] (linear)
        a_value = 0.4;
        T = @(x,n) a_value.^(n-1) .* x .* (n >= 0);
    elseif s == 3
        % System 3: y[n] = x[n] + x^2[n] (non-linear)
        T = @(x,n) x + x.^2;
    elseif s == 4
        % System 4: y[n] = x[n] + 10 * u[n-1] (non-linear)
        T = @(x,n) x + 10 * (n >= 1);
    end
    
    % Test linearity
    a = 5; % Scaling factor
    b = -0.5;

    y1 = T(x1, n); % Output for x1
    y2 = T(x2, n); % Output for x2

    % Check if T{a*x1 + b*x2} equals a*y1 + b*y2 (linearity condition)
    LHS = T(a * x1 + b * x2, n); % Left-hand side
    RHS = a * y1 + b * y2; % Right-hand side

    % Plot results to compare
    figure;
    subplot(2,1,1);
    stem(n, LHS, 'b', 'filled');
    title('LHS: T\{a*x1 + b*x2\}');
    xlabel('n');
    ylabel('Amplitude');
    
    subplot(2,1,2);
    stem(n, RHS, 'r', 'filled');
    title('RHS: a*y1 + b*y2');
    xlabel('n');
    ylabel('Amplitude');
  
    % Pause for inspection
    disp('Paused for inspection...');
    pause;
end

%% Test for Time Invariance

% Define the unit step function u[n]
u = (n >= 0); % unit step function (u[n] = 1 for n >= 0, 0 otherwise)

% Define the input signal in the form of a^n * u[n]
a1 = 0.9;
x = a1.^n .* u; % Input signal

% Define time shift value
n_shift = 3; % Shift the signal by 3 units

% Shifted input
x_shifted = [zeros(1, n_shift), x(1:end-n_shift)]; % Right shift by n_shift

% Define the systems
for s = 1:1:4
    if s == 1
        % System 1: y[n] = x[n] + x[n+5] (time invariant)
        T = @(x,n) x + [x(6:end) zeros(1,5)];
    elseif s == 2
        % System 2: y[n] = a^(n-1) * x[n] (time varying)
        a_value = 0.4;
        T = @(x,n) a_value.^(n-1) .* x .* (n >= 0);
    elseif s == 3
        % System 3: y[n] = x[n] + x^2[n] (time invariant)
        T = @(x,n) x + x.^2;
    elseif s == 4
        % System 4: y[n] = x[n] + 5 * u[n-1] (time varying)
        T = @(x,n) x + 5 * (n >= 1);
    end

    % Original system output
    y_original = T(x, n);

    % Output for shifted input
    y_shifted_input = T(x_shifted, n);

    % Shift the output of the original input
    y_original_shifted = [zeros(1, n_shift), y_original(1:end-n_shift)];

    % Plot results to compare
    figure;
    subplot(2,1,1);
    stem(n, y_shifted_input, 'b', 'filled');
    title('T{[x[n-n\_shift]]} (System applied to shifted input)');
    xlabel('n');
    ylabel('Amplitude');

    subplot(2,1,2);
    stem(n, y_original_shifted, 'r', 'filled');
    title('y[n-n\_shift] (Shifted output of original input)');
    xlabel('n');
    ylabel('Amplitude');

    disp('Paused for inspection...');
    pause;
end

%% Impulse Response Derivation

% Define a system (Example: Moving average system)
T = @(x) filter([1, 1, 1], 1, x); % 3-point moving average

% Define a unit impulse signal
impulse = (n == 0); % Unit impulse delta[n]

% Derive impulse response
h = T(impulse); % Impulse response is the output for the impulse input

% Plot impulse response
figure;
stem(n, h, 'filled');
title('Impulse response of the LTI system');
xlabel('n');
ylabel('Amplitude');

disp('Paused for inspection...');
pause;

% Define another system
T = @(x) x + [x(6:end), zeros(1, 5)];

% Derive impulse response
h = T(impulse); % Impulse response is the output for the impulse input

% Plot impulse response
figure;
stem(n, h, 'filled');
title('Impulse response of the LTI system');
xlabel('n');
ylabel('Amplitude');

disp('Paused for inspection...');
pause;

%% Convolution

% Define the unit step function
u = (n >= 0); % Unit step function (u[n] = 1 for n >= 0, 0 otherwise)

% Define the input signal: a^n * u[n]
a1 = 0.9; % Base for exponential input
x = a1.^n .* u; % Input signal: a1^n * u[n]

% Define system: y[n] = x[n] + x[n+5] (LTI)
h = [1 zeros(1, 4) 1]; % Impulse response

% Perform convolution
y_conv = conv(x, h);

% Truncate the convolution result to match the input signal length
y_conv_truncated = y_conv(1:length(n));

% Plot input and truncated convolution output
figure;
subplot(2,1,1);
stem(n, x, 'b', 'filled');
title('Input signal x[n]');
xlabel('n');
ylabel('Amplitude');

subplot(2,1,2);
stem(n, y_conv_truncated, 'r', 'filled');
title('Output y[n] = x[n] + x[n+5] (truncated convolution result)');
xlabel('n');
ylabel('Amplitude');s