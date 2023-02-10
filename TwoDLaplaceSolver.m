% Parameters for solution var
nx = 100;
ny = 100;
maxNumberInter = 10000;

% V Matrix as a solution vairable
V = zeros(nx, ny);

% Initial BC for the solution
V(:, 1) = 1;
V(:, ny) = 0;
V(1, :) = V(2, :);
V(ny, :) = V(ny-1, :);

% Interation startes from here
for temp = 1:maxNumberInter
    % Get new solution
    for i = 2:nx-1
        for j = 2:ny-1
            V(i, j) = (V(i-1, j) + V(i+1, j) + V(i, j-1) + V(i, j+1))/4;
        end
    end
    
    % Reset boundary conditions

    V(1, :) = V(2, :);
    V(ny, :) = V(ny-1, :);
end

% Plot final solution using surf()
surf(V);

% Calculates electric field
[Ex, Ey] = gradient(-V);

% Electric field plot using surf() and quiver()
figure;
surf(Ex);
figure;
surf(Ey);
figure;
quiver(Ex, Ey);

% Image processing using imboxfilt()
%V_filt = imboxfilt(V, 3);

%filtered solution Plot 
%figure;
%surf(V_filt);