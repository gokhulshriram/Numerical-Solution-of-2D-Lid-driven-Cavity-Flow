clear;
close all;

%% Variables
Re = 0.1;
uwall = 1;
max_time = 0.1;
dt = 0.00001;
L = 1;
h = 1/50;
x = 0:h:L;
y = 0:h:L;
timecounter = 0;
nx = L/h +1;
ny = nx;
u = zeros(nx,nx);
v = u;
omega = u;
psi = u;
p = u;
p_old = p;
[X, Y] = meshgrid(x, y);

%% Fluid Flow Simulation

while (timecounter < max_time)
    % Set boundary conditions
    for i = 1:nx
        omega(i, end) = -2 * psi(i, end - 1) / (h * h) - uwall * 2 / h;
        omega(i, 1) = -2 * psi(i, 2) / (h * h);
        omega(1, i) = -2 * psi(2, i) / (h * h);
        omega(end, i) = -2 * psi(end - 1, i) / (h * h);
    end

    omega_old = omega;
    psi_old = psi;
    p_old = p;

    for i = 2:nx - 1
        for j = 2:ny - 1
            omega(i, j) = omega_old(i, j) + ((psi_old(i, j - 1) - psi_old(i, j + 1)) * (omega_old(i - 1, j) - omega_old(i + 1, j)) / (4 * h * h) + (psi_old(i + 1, j) - psi_old(i - 1, j)) * (omega_old(i, j + 1) - omega_old(i, j - 1)) / (4 * h * h) + (1 / Re) * (omega_old(i + 1, j) + omega_old(i, j + 1) - 4 * omega_old(i, j) + omega_old(i - 1, j) + omega_old(i, j - 1)) / (h * h)) * dt;
            psi(i, j) = (1 / 4) * (omega_old(i, j) * h * h + psi_old(i + 1, j) + psi_old(i, j + 1) + psi_old(i, j - 1) + psi_old(i - 1, j));
            u(i, j) = (psi(i, j + 1) - psi(i, j - 1)) / (2 * h);
            v(i, j) = -(psi(i + 1, j) - psi(i - 1, j)) / (2 * h);

            %p(i,j) = (h^2 / 4) * ( (p_old(i+1,j) + p_old(i-1,j) + p_old(i,j+1) + p_old(i,j-1))/(h^2) + ( (u(i+1,j) - u(i-1,j))/(2*h) )^2 + ( v(i,j+1) - v(i,j-1) * (1/(2*h)) )^2 + 2 * ( (u(i,j+1) - u(i,j-1)) * (v(i+1,j) - v(i-1,j)) ) / (4*h^2) );
        end
    end
    %     if rem(timecounter,0.02) == 0
    %         figure;
    %         streamline(X,Y,u',v');
    %         xlabel('X');
    %         ylabel('Y');
    %         title('Streamline at t=0.02')
    %         axis equal;
    %     end
    timecounter = timecounter + dt;
end
u(2:nx,end) = uwall;
u1 = u;
v1 = v;

%for Re = 1
Re = 1;
dt = 0.0001;

while (timecounter < max_time)
    % Set boundary conditions
    for i = 1:nx
        omega(i, end) = -2 * psi(i, end - 1) / (h * h) - uwall * 2 / h;
        omega(i, 1) = -2 * psi(i, 2) / (h * h);
        omega(1, i) = -2 * psi(2, i) / (h * h);
        omega(end, i) = -2 * psi(end - 1, i) / (h * h);
    end

    omega_old = omega;
    psi_old = psi;
    p_old = p;

    for i = 2:nx - 1
        for j = 2:ny - 1
            omega(i, j) = omega_old(i, j) + ((psi_old(i, j - 1) - psi_old(i, j + 1)) * (omega_old(i - 1, j) - omega_old(i + 1, j)) / (4 * h * h) + (psi_old(i + 1, j) - psi_old(i - 1, j)) * (omega_old(i, j + 1) - omega_old(i, j - 1)) / (4 * h * h) + (1 / Re) * (omega_old(i + 1, j) + omega_old(i, j + 1) - 4 * omega_old(i, j) + omega_old(i - 1, j) + omega_old(i, j - 1)) / (h * h)) * dt;
            psi(i, j) = (1 / 4) * (omega_old(i, j) * h * h + psi_old(i + 1, j) + psi_old(i, j + 1) + psi_old(i, j - 1) + psi_old(i - 1, j));
            u(i, j) = (psi(i, j + 1) - psi(i, j - 1)) / (2 * h);
            v(i, j) = -(psi(i + 1, j) - psi(i - 1, j)) / (2 * h);

            %p(i,j) = (h^2 / 4) * ( (p_old(i+1,j) + p_old(i-1,j) + p_old(i,j+1) + p_old(i,j-1))/(h^2) + ( (u(i+1,j) - u(i-1,j))/(2*h) )^2 + ( v(i,j+1) - v(i,j-1) * (1/(2*h)) )^2 + 2 * ( (u(i,j+1) - u(i,j-1)) * (v(i+1,j) - v(i-1,j)) ) / (4*h^2) );
        end
    end
    %     if rem(timecounter,0.02) == 0
    %         figure;
    %         streamline(X,Y,u',v');
    %         xlabel('X');
    %         ylabel('Y');
    %         title('Streamline at t=0.02')
    %         axis equal;
    %     end
    timecounter = timecounter + dt;
end
%% Calculate pressure
u(2:nx,end) = uwall;
for i = 2:nx - 1
    for j = 2:ny - 1
        p(i,j) = (h^2 / 4) * ( (p_old(i+1,j) + p_old(i-1,j) + p_old(i,j+1) + p_old(i,j-1))/(h^2) + ( (u(i+1,j) - u(i-1,j))/(2*h) )^2 + ( v(i,j+1) - v(i,j-1) * (1/(2*h)) )^2 + 2 * ( (u(i,j+1) - u(i,j-1)) * (v(i+1,j) - v(i-1,j)) ) / (4*h^2) );
    end
end

%% Plot velocity contours
% figure;
% subplot(1,2,1);
% contourf(X, Y, u' , 20, 'LineStyle', 'none');
% colorbar
% title('Velocity Contours');
% xlabel('x');
% ylabel('y');
% colormap('turbo');
% axis equal;
% subplot(1,2,2);
% contourf(X , Y , v' , 20, 'LineStyle', 'none');
% colorbar;
% title('v Velocity Contours');
% xlabel('x');
% ylabel('y');
% colormap('turbo');
% axis equal;
% 
% figure;
% subplot(1,2,1);
% streamslice(X,Y,u',v');
% title('Streamlines');
% xlabel('X');
% ylabel('Y');
% xlim([0 1]);
% ylim([0 1]);
% axis equal;
% axis tight;
% subplot(1,2,2);
% contourf(X, Y, p' , 20, 'LineStyle', 'none');
% colorbar
% title('Pressure');
% xlabel('x');
% ylabel('y');
% colormap('turbo');
% axis equal;
% 


u = u';
v = v';
u1 = u1';
v1 = v1';

figure;
subplot(1,2,1);
plot(y,u1(:,round(nx/2)),'bo-', 'LineWidth', 2);
hold on;
plot(y,u(:,round(nx/2)),'r--', 'LineWidth', 2);
legend('Re = 0.1','Re = 1');
xlabel('y');
ylabel('Velocity');
title('U velocity');
axis tight;
hold off;

subplot(1,2,2);
plot(y,v1(:,round(nx/2)),'bo-', 'LineWidth', 2);
hold on;
plot(y,v(:,round(nx/2)),'r--', 'LineWidth', 2);
legend('Re = 0.1','Re = 1');
xlabel('y');
ylabel('Velocity');
title('V velocity');
axis tight;
hold off;

% function p = calcPress(u,v,h,p_old)
% u(2:nx,end) = uwall;
% for i = 2:nx - 1
%     for j = 2:ny - 1
%         p(i,j) = (h^2 / 4) * ( (p_old(i+1,j) + p_old(i-1,j) + p_old(i,j+1) + p_old(i,j-1))/(h^2) + ( (u(i+1,j) - u(i-1,j))/(2*h) )^2 + ( v(i,j+1) - v(i,j-1) * (1/(2*h)) )^2 + 2 * ( (u(i,j+1) - u(i,j-1)) * (v(i+1,j) - v(i-1,j)) ) / (4*h^2) );
%     end
% end
% end
