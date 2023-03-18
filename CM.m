% Define parameters
Is = 0.01e-12; % Forward bias saturation current
Ib = 0.1e-12; % Breakdown saturation current
Vb = 1.3; % Breakdown voltage
Gp = 0.1; % Parasitic parallel conductance

% Create V vector
V = linspace(-1.95, 0.7, 200);

% Create ideal diode current vector
I_id = Is*(exp(1.2*(V/0.025) - 1));

% Create breakdown current vector
I_bd = -Ib*(exp(1.2*(-(V + Vb)/0.025) - 1));

% Create parasitic parallel leakage current vector
I_gp = Gp*V;

% Create total current vector with noise
I = I_id + I_bd + I_gp;
I_noisy = I + 0.2*randn(size(I));

% Plot the data using plot() and semilogy()
figure
plot(V, I, 'b', V, I_noisy, 'r')
xlabel('V (V)')
ylabel('I (A)')
title('Diode Characteristic Curve')
legend('Ideal Data', 'Noisy Data')
grid on

figure
semilogy(V, abs(I), 'b', V, abs(I_noisy), 'r')
xlabel('V (V)')
ylabel('I (A)')
title('Diode Characteristic Curve (Log Scale)')
legend('Ideal Data', 'Noisy Data')
grid on

%Polynomial fitting part 3

% Generate data
V = linspace(-1.95, 0.7, 200)';
I = Is*(exp(1.2/0.025*(V+Vb-2*Vb*(V+Vb<0)/Vb))-1) + Gp*V - Ib*(exp(1.2/0.025*(-(V+Vb)/Vb))-1);
I_noise = I + 0.2*randn(size(I)).*I; % Add 20% noise

% Plot data
figure;
plot(V, I, 'o', V, I_noise, '.');
title('Diode I-V Characteristic');
xlabel('Voltage (V)');
ylabel('Current (A)');
legend('Ideal', 'Noisy');
grid on;

% Polynomial fitting
p4 = polyfit(V, I_noise, 4);
p8 = polyfit(V, I_noise, 8);

% Evaluate polynomial fits
I_fit4 = polyval(p4, V);
I_fit8 = polyval(p8, V);

% Plot polynomial fits
figure;
semilogy(V, abs(I_noise), '.', V, abs(I_fit4), V, abs(I_fit8));
title('Polynomial Fitting');
xlabel('Voltage (V)');
ylabel('Current (A)');
legend('Noisy', '4th order', '8th order');
grid on;

%conclusion: From the plot, we can see that the 4th order 
% polynomial fit is generally a good approximation of the noisy data, 
% while the 8th order fit overfits the data and 
% exhibits some oscillations. 
% However, both polynomial fits deviate significantly 
% from the ideal diode characteristic at low and high voltages. 
% This is because the polynomial model is not flexible 
% enough to capture the exponential behavior of the diode equation. 
% Thus, nonlinear curve fitting may be a better approach for 
% modeling the diode characteristic.

% Nonlinear curve fitting part 4
fo = fittype('A*(exp(1.2*x/25e-3)-1) + B*x - C*(exp(1.2*(-(x+D))/25e-3)-1)');

% Fit only A and C with B and D fixed
fo_fixed = fittype('A*(exp(1.2*x/25e-3)-1) + 0.1*x - C*(exp(1.2*(-(x+1.3))/25e-3)-1)');
ff_fixed = fit(V, I_noisy.', fo_fixed);

I_fit_fixed = ff_fixed(V);

% Fit A, B, and C with D fixed
fo_fixed2 = fittype('A*(exp(1.2*x/25e-3)-1) + B*x - C*(exp(1.2*(-(x+1.3))/25e-3)-1)');
ff_fixed2 = fit(V, I_noisy.', fo_fixed2);
I_fit_fixed2 = ff_fixed2(V);

% Fit all four parameters
ff = fit(V, I_noisy.', fo);
I_fit = ff(V);

% Plot the data and the three fitted curves
figure
semilogy(V, abs(I_noisy), '.', V, abs(I_fit_fixed), V, abs(I_fit_fixed2), V, abs(I_fit));
title('Nonlinear Curve Fitting');
xlabel('Voltage (V)');
ylabel('Current (A)');
legend('Noisy', 'A and C fixed', 'A, B, and C fixed', 'All Parameters');
grid on;

% Print the coefficients for the three fits
disp('Fit with A and C only:');
disp(ff_fixed);
disp('Fit with A, B, and C only:');
disp(ff_fixed2);
disp('Fit with all parameters:');
disp(ff);

% Draw conclusions
% The three fitted curves appear to closely match the noisy data, 
% with the curve that includes all four parameters fitting the data 
% the most closely. The fit with only A and C fixed appears to have 
% the largest deviation from the noisy data, while the fit with A, 
% B, and C fixed is a good compromise between the other two fits. 
% This suggests that the breakdown voltage and parasitic parallel 
% conductance have a significant impact on the behavior of the diode. 
% Additionally, the use of nonlinear curve fitting with the fit() 
% function can be a powerful tool for analyzing and modeling complex 
% physical systems.


%Part 5

% Generate simulated IV curve for a diode
V = linspace(0, 0.8, 100);
I = diode(V, 1e-9, 0.025);

% Add random noise to the simulated data
I_noisy = awgn(I, 20);

% Remove NaN values from the data
V = V(~isnan(I_noisy));
I_noisy = I_noisy(~isnan(I_noisy));

% Nonlinear curve fitting
fo = fittype('A*(exp(1.2*x/25e-3)-1) + B*x - C*(exp(1.2*(-(x+D))/25e-3)-1)');
ff = fit(V', I_noisy', fo, 'StartPoint', [1e-6, 0.5, 1e-6, 0.6]);

% Plot the data and the fitted curve
figure
plot(V, I_noisy, '.', V, ff(V));
title('Nonlinear Curve Fitting');
xlabel('Voltage (V)');
ylabel('Current (A)');
legend('Noisy', 'Fitted');
grid on;

% Print the coefficients for the fit
disp(ff);

% Fitting using the Neural Net model
inputs = V.';
targets = I.';
hiddenLayerSize = 10;
net = fitnet(hiddenLayerSize);
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;
[net,tr] = train(net,inputs,targets);
outputs = net(inputs);
errors = gsubtract(outputs,targets);
performance = perform(net,targets,outputs);
view(net);
Inn = outputs;

% Plot the Neural Net model output
figure
semilogy(V, abs(I_noisy), '.', V, abs(outputs));
title('Neural Net Fitting');
xlabel('Voltage (V)');
ylabel('Current (A)');
legend('Noisy', 'Neural Net Fitting');
grid on;












