clear all; 
close all;

a = 1; %side a is 1m long
V0 = 1; %rectangle is at potential of 1V
EPS0 = 8.8542*10^(-12);
n = 10; %number of patches per side a
m = 20; %added m = 20 

%coordinates of centers of patches

dp = a/n; %size of individual patch

for i = 1:n
    for j = 1:m %M replaced N in this for loop
        x(m*(i-1)+j) = (i-1/2)*dp;
        y(m*(i-1)+j) = (j-1/2)*dp;
        dS(m*(i-1)+j) = dp^2;
    end
end

%matrix A 
N = length(x);
for i = 1 : N
    for j = 1 : N
        r = sqrt ((x(j)-x(i))^2 + (y(j)-y(i))^2);
        if (i==j)
            A(i,j) = sqrt(dS(j))/(2*sqrt(pi)*EPS0);
        else
            A(i,j) = dS(j)/(4*pi*EPS0*r);   
        end
    end
end

%matrix B

B = V0*ones(1,N)';  %transpose -- to make it a column matrix

rhos = A\B;  %solving the equation

% making 2D matrix rhos2D from results rhos in 1D array

for i = 1:N
    row_index = ceil(x(i)/dp);
    column_index = ceil(y(i)/dp);
    rhos2D(row_index, column_index) = rhos(i);
end

figure(1);

[x2D,y2D] = meshgrid(dp/2:dp:a, dp/2:dp:2*a);%Changed to 2*a for y-axis
surf(x2D, y2D, rhos2D');
colormap('cool');
shading interp;
title(' Surface charge distribution of the 1 X 2 rectangle plate');
xlabel('x [m]');
ylabel('y [m]');
zlabel('\rho_s [^{C}/_{m^2}]');


figure(2);

pcolor(x2D,y2D,rhos2D');
colormap('cool');
shading interp;
title('Surface charge density of the 1 X 2 rectangle plate in C/m^2');
xlabel('x');
ylabel('y');

Qtot = dS*rhos;

disp(['Total charge of  is 1 X 2 rectangle is ' num2str(Qtot*10^12) ' pC']);