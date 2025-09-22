clear all; 
close all;

a = 1; %side a is 1m long
V0 = 1; %rectangle is at potential of 1V
EPS0 = 8.8542*10^(-12);
n = 10; %number of patches per side a

%coordinates of centers of patches

dp = a/n; %size of individual patch
x = [];
y = [];
ds = [];
for i = 1:(2*n)
    for j = (n+1):(2*n)
        yj = (j-1)*dp; % y-coordinate of the patch center
        xi = (i-1)*dp; % x-coordinate of the patch center
        x(end+1) = xi;
        y(end+1) = yj;
        ds(end+1) = dp^2;
    end
end

for i = (n+1):(3*n)
    for j = 1:n
        yj = (j-1)*dp; % y-coordinate of the patch center
        xi = (i-1)*dp; % x-coordinate of the patch center
        x(end+1) = xi;
        y(end+1) = yj;
        ds(end+1) = dp^2;
    end
end

%matrix A 
N = length(x); % Calculate the number of patches based on x coordinates
A = zeros(N,N);
for i = 1 : N
    for j = 1 : N
        if (i==j)
            A(i,j) = sqrt(ds(j))/(2*sqrt(pi)*EPS0);
        else
            rij = hypot(x(j)-x(i), y(j)-y(i));
            A(i,j) = ds(j)/(4*pi*EPS0*rij);   
        end
    end
end

%matrix B

B = V0*ones(1,N)';  %transpose -- to make it a column matrix

rhos = A\B;  %solving the equation
rhos2D = NaN(3*n,2*n);

% making 2D matrix rhos2D from results rhos in 1D array

for k = 1:N
    row_index = floor(x(k)/dp) + 1;
    column_index = floor(y(k)/dp) + 1;
    rhos2D(row_index, column_index) = rhos(k);
end

figure(1);

[x2D,y2D] = meshgrid(0:dp:3*a-dp,0:dp:2*a-dp);
surf(x2D, y2D, rhos2D');
colormap('cool');
shading interp;
title(' Surface charge distribution of 3 X 2 rectangle');
xlabel('x [m]');
ylabel('y [m]');
zlabel('\rho_s [^{C}/_{m^2}]');


figure(2);

pcolor(x2D,y2D,rhos2D');
colormap('cool');
shading interp;
title('Surface charge density of 3 X 2 rectangle in C/m^2');
xlabel('x');
ylabel('y');

Qtot = ds*rhos;

disp(['Total charge of 3 x 2 rectangle is ' num2str(Qtot*10^12) ' pC']);