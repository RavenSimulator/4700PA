% PA7
% Parameters
R1 = 1; G1 = 1/R1;
C1 = 0.25;
R2 = 2; G2 = 1/R2;
L = 0.2;
R3 = 10; G3 = 1/R3;
alpha = 100;
R4 = 0.1; G4 = 1/R4;
R0 = 1000; G0 = 1/R0;
w = linspace(0,100);

%F = [V1; V2; V3; V4; V0; Is; Is2; IL]

% G matrix
G = [1 0 0 0 0 0 0 0;
    G1 -G1 0 0 0 -1 0 0 ;
    G1 -G1-G2 0 0 0 0 0 -1;
    0 0 -G3 0 0 0 0 1;
    0 1 -1 0 0 0 0 0;
    0 0 0 G4 -G4 0 1 0;
    0 0  -alpha*G3 1 0 0 0 0;
    0 0 0 -G4 G4-G0 0 0 0];

% C matrix
C = [0 0 0 0 0 0 0 0 ;
     C1 -C1 0 0 0 0 0 0;          
     C1 -C1 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 ;
     0 0 0 0 0 0 0 L;
     0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0];

%A = G + 1i*x*C;

%F = [Vin; 0; 0; 0; 0; 0; 0; 0];
Fw = [1; 0; 0; 0; 0; 0; 0; 0];

% Plot 1
i = 1;
for Vin = -10:10
    F = [Vin; 0; 0; 0; 0; 0; 0; 0];
    x = inv(G)*F;

    y_v3(i) = x(3);
    y_v5(i) = x(5);

    i = i+1;

end

% Plot 2
j = 1;
for w = 0:100
    Vin = 1;

    F = [Vin; 0; 0; 0; 0; 0; 0; 0];
    A = G + 1i*w*C;

    x = inv(A)*F;

    gain_3(j) = abs(x(3)/Vin);
    gain_0(j) = abs(x(5)/Vin);

    j = j + 1;
end

% MC
C_rand = 0.05*randn(1000)+0.25;

for i = 1:size(C_rand)
    Vin = 1;
    w = 3*pi;
    F = [Vin; 0; 0; 0; 0; 0; 0; 0];

    C_ran = [0 0 0 0 0 0 0 0 ;
             C_rand(i) -C_rand(i) 0 0 0 0 0 0;
             C_rand(i) -C_rand(i) 0 0 0 0 0 0;
             0 0 0 0 0 0 0 0;
             0 0 0 0 0 0 0 0;
             0 0 0 0 0 0 0 0;
             0 0 0 0 0 0 0 0;
             0 0 0 0 0 0 0 0];

    A = G + 1i*w*C_ran;
    x = inv(A)*F;

    gain_C(i) = abs(x(5)/Vin);
end

Vin = -10:10;
w = 0:100;

figure()
subplot(2,2,1)
plot(Vin,y_v3);hold on
plot(Vin,y_v5);hold off
xlabel('Vin')
ylabel('V')
legend('V3','V_{out}')

subplot(2,2,2)
plot(w,gain_3);hold on
plot(w,gain_0);hold off
xlabel('w')
ylabel('Gain')
legend('V3','V_{out}')

subplot(2,2,3)
histogram(C_rand);
xlabel('C')
ylabel('#')

subplot(2,2,4)
histogram(gain_C);
xlabel('gain')
ylabel('#')