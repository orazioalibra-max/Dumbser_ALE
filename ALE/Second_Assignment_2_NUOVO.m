clc; clear; close all;

%% Parameters
gamma = 1.4; % Ratio of specific heats
Nx = 500;    % Number of points
x_min = -1; x_max = 1; % Domain limits
CFL = 0.9;   % CFL
t_final = 0.4; % Final time

% Discretization
dx = (x_max - x_min) / Nx;
xl = linspace(x_min, x_max, Nx+1); % Interface positions
x = (xl(1:end-1) + xl(2:end)) / 2; % Cell centers

% Initial conditions
rhoL = 1.0; uL = 0.0; pL = 1.0;  % Left state
rhoR = 0.125; uR = 0.0; pR = 0.1; % Right state

% Compute initial conserved variables
E_L = pL / (gamma - 1) + 0.5 * rhoL * uL^2;
E_R = pR / (gamma - 1) + 0.5 * rhoR * uR^2;

% Initialize Q (conserved variables)
N_half = floor(Nx/2); % removes division problems
Q = zeros(3, Nx);
Q(:, 1:N_half) = repmat([rhoL; rhoL * uL; E_L], 1, N_half);
Q(:, N_half+1:Nx) = repmat([rhoR; rhoR * uR; E_R], 1, Nx - N_half);

%% Time stepping
t = 0;
while t < t_final

    % Compute how large are the cell
    dx = xl(2:end) - xl(1:end-1);

    % Compute velocity of moving nodes
    vel = zeros(1, Nx+1);
    for i = 2:Nx
        u_i = Q(2, i) / Q(1, i);
        u_im1 = Q(2, i-1) / Q(1, i-1);
        vel(i) = (dx(i) * u_i + dx(i-1) * u_im1 ) / (dx(i-1) + dx(i));
    end

    % Boundary values (transmissive)
    vel(1) = vel(2);
    vel(Nx+1) = vel(Nx);

    % Compute the time step using CFL condition
    dt=1e-2;
    for i=1:Nx
        a_max = max( max( abs(lambda(Q(:,i),gamma)-vel(i+1))), max( abs(lambda(Q(:,i),gamma)-vel(i))));
        dt = min(dt, CFL*dx(i)/a_max);
    end 

    if t + dt > t_final
        dt = t_final - t;
    end

    % Compute new interface positions
    xl_new = xl + dt*vel; %for stability

    % Compute new cell centers and dx
    x_new = (xl_new(1:end-1) + xl_new(2:end)) / 2;
    dx_new = xl_new(2:end) - xl_new(1:end-1);

    % Compute fluxes and max wave speed
    Q_new=Q;    
    for i = 2:Nx-1
        Q_L = Q(:, i-1);
        Q_i = Q(:, i);
        Q_R = Q(:, i+1);
        
        % length of the segment of the control space-time volume
        a  = dx(i);
        c  = dx_new(i);
        b  = dt*sqrt(1+vel(i+1)^2);
        d  = dt*sqrt(1+vel(i)^2);
        n_R = [1,-vel(i+1)]/(sqrt(1+vel(i+1)^2));
        n_L = [-1,vel(i)]/(sqrt(1+vel(i)^2));

        % Compute ALE Rusanov flux
        smaxL = max( max( abs(lambda(Q_L,gamma)-vel(i))), max( abs(lambda(Q_i,gamma)-vel(i))));
        smaxR = max( max( abs(lambda(Q_i,gamma)-vel(i+1))), max( abs(lambda(Q_R,gamma)-vel(i+1))));
        fiL = 0.5*(F(Q_i,gamma) + F(Q_L,gamma))*(n_L(1)) + 0.5*(Q_i+Q_L)*(n_L(2)) - 0.5*smaxL*(Q_i-Q_L)*(n_L(1));
        fiR = 0.5*(F(Q_i,gamma) + F(Q_R,gamma))*(n_R(1)) + 0.5*(Q_R+Q_i)*n_R(2) - 0.5*smaxR*(Q_R-Q_i)*(n_R(1));
        Q_new(:,i) = dx(i)/dx_new(i)*Q_new(:,i) - 1/dx_new(i)*(b*fiR+d*fiL);
    end
    
    % Apply transmissive boundary conditions
    Q_new(:, 1) = Q_new(:, 2);
    Q_new(:, Nx) = Q_new(:, Nx-1);
    
    % Update variables
    Q = Q_new;
    xl = xl_new;
    x = x_new;
    dx = dx_new;
    t = t + dt;

    % interpolate data from the cells to the nodes
    TotArea = zeros(Nx+1,1); % Total area of the cells attached to a node
    QNode   = zeros(3,Nx+1); % Averaged state in the node
    for i = 1:Nx
        for k = 1:2
            iNode = i+k-1; % Global node number
            TotArea(iNode) = TotArea(iNode) + dx_new(i);
            QNode(:,iNode) = QNode(:,iNode) + dx_new(i)*Q(:,i);
        end
    end

    % area-weighted average
    for i = 1:Nx+1
        QNode(:,i) = QNode(:,i)/TotArea(i);
    end

    %Plot results
    rho = QNode(1, :);
    u = QNode(2, :) ./ rho;
    E = QNode(3, :);
    p = (gamma - 1) * (E - 0.5 * rho .* u.^2);
    
    subplot(2,2,1); plot(xl, rho, 'bo'); title(['Density at t = ', num2str(t)]); xlabel('x'); ylabel('\rho');
    subplot(2,2,2); plot(xl, u, 'ro'); title('Velocity'); xlabel('x'); ylabel('u');
    subplot(2,2,3); plot(xl, p, 'go'); title('Pressure'); xlabel('x'); ylabel('p');
    subplot(2,2,4); plot(xl, E, 'mo'); title('Energy'); xlabel('x'); ylabel('E');
    drawnow;
end

%% Function
function flux = F(Q,gamma) 
% Flux tensor 

Q = Q(:);
u = Q(2) / Q(1);

p = (gamma-1) * ( Q(3) - 0.5*Q(1)*(u^2) );

flux = [Q(2); 
        u * Q(2) + p;
        u * (Q(3) + p)];
end

%% Function: Compute eugenvalues of the PDE
function L=lambda(Q,gamma)

u = Q(2)/Q(1); 

p = (gamma-1) * ( Q(3) - 0.5*Q(1)*(u^2) );
c = sqrt(gamma*p/Q(1));

L = [u-c, u, u+c];
end