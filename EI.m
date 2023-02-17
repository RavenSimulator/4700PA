% Set the grid size
nx = 50;
ny = 50;

% Set the grid spacing
L = 1;
dx = L/(nx-1);
dy = L/(ny-1);

% Set the time step
dt = 0.001;

% Create the G matrix
G = sparse(nx*ny, nx*ny);

% Set the boundary nodes
for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;
        if i == 1 || i == nx || j == 1 || j == ny
            G(n,n) = 1;

        elseif(i>10 & i<20 & j>10 & j<20)
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;

            G(n,n) = -2;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nym) = 1;
            G(n,nyp) = 1; 
        
        else
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;

            G(n,n) = -4;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nym) = 1;
            G(n,nyp) = 1;
            
        end 
    end
end

% Plottig the G matrix
figure;
spy(G);

% Compute the eigenvectors and eigenvalues
num_eig = 9;
options.issym = true;
options.isreal = true;
[E,D] = eigs(G, num_eig, 'SM', options);

% Plot the eigenvalues
figure;
plot(diag(D), 'o');
title('Eigenvalues');

% Plot the eigenvectors
figure;
for i = 1:num_eig
    subplot(3, 3, i);
    eig_vec = E(:,i);
    eig_mat = reshape(eig_vec, [ny, nx]);
    surf(eig_mat');
    title(['Eigenvector ', num2str(i)]);
    xlabel('x');
    ylabel('y');
end
