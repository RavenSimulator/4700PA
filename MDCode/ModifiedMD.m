% Parameters
n_atoms = 10;        % Number of atoms
k = 1.0;             % Spring constant
x0 = 1.0;            % Equilibrium position of spring
dt = 0.01;           % Time step
t_max = 10.0;        % Maximum simulation time
sigma = 0.1;         % Initial perturbation amplitude
box_size = 10.0;     % Size of simulation box
atom_radius = 0.1;   % Radius of atoms

% Initialize positions and velocities
x = rand(1, n_atoms) * box_size;
v = zeros(size(x));
x = x + sigma * randn(size(x));

% Main loop
for t = 0:dt:t_max
    % Calculate forces
    f = zeros(size(x));
    for i = 2:n_atoms-1
        f(i) = k * ((x(i+1) - x(i)) + (x(i-1) - x(i)) - 2*x0);
    end
    % Apply boundary conditions (periodic)
    f(1) = k * ((x(2) - x(1)) + (x(n_atoms) - x(1)) - 2*x0);
    f(n_atoms) = k * ((x(1) - x(n_atoms)) + (x(n_atoms-1) - x(n_atoms)) - 2*x0);
    % Update positions and velocities
    x_new = x + v*dt + 0.5*f*dt^2;
    v_new = v + 0.5*f*dt;

    % Apply boundary conditions (bounce)
    for i = 1:n_atoms
        if x_new(i) < atom_radius
            v_new(i) = abs(v_new(i));
            x_new(i) = atom_radius;
        elseif x_new(i) > box_size - atom_radius
            v_new(i) = -abs(v_new(i));
            x_new(i) = box_size - atom_radius;
        end
        % Check for collisions with other atoms
        for j = i+1:n_atoms
            if abs(x_new(i) - x_new(j)) < 2*atom_radius
                v_temp = v_new(i);
                v_new(i) = v_new(j);
                v_new(j) = v_temp;
                % Move atoms apart if they overlap
                if x_new(i) < x_new(j)
                    x_new(i) = x_new(j) - 2*atom_radius;
                else
                    x_new(j) = x_new(i) - 2*atom_radius;
                end
            end
        end
    end
    x = x_new;
    v = v_new + 0.5*f*dt;
    % Plot positions
    plot(x, zeros(size(x)), 'o');
    xlim([0, box_size]);
    pause(0.1);
end
