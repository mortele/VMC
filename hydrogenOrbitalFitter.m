%% Clean up
close all;
clear variables;
clc;
format;


%% 
N = 1000;
alpha = 1; %3.983;
x = linspace(0,3,N);
Z = 1;
n = 1;
l = 0;
m = 0;

norm = sqrt(4*alpha.^3); %((2*Z)/n)^(1/4) * sqrt(factorial(n-l-1)/(2*n*factorial(n+l)));

f = @(x) norm .* exp(-alpha.*x);
y = f(x);



nor = @(a) sqrt(4*pi)*(2*a/pi).^(3./4);
G   = @(c,a,x) nor(a) * c * exp(-a.*x.^2); 

g = @(x) ( G(0.15432897,30.1678710,x) + G(0.53532814,5.4951153,x) + G(0.44463454,1.4871927,x));
h = @(x) (G(0.00916359628,312.87049370,x) + G(0.04936149294,57.364462530,x) + G(0.16853830490,16.048509400,x) + G(0.37056279970,5.5130961190,x) + G(0.41649152980,2.1408965530,x) + G(0.13033408410,0.8817394283,x));

dx = mean(diff(x));
If = sum(dx * x.^2.*f(x).^2)hydrog
Ig = sum(dx * x.^2.*g(x).^2)
Ih = sum(dx * x.^2.*h(x).^2)

plot(x,x.^2.*f(x),'k--');
grid('on')
hold('on')
plot(x,x.^2.*g(x)./Ig,'b-')
plot(x,x.^2.*h(x)./Ih,'r-')
h = legend('Hydrogen','STO-3G','STO-6G');
set(h, 'FontSize', 16,'interpreter','latex');