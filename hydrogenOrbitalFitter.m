%% Clean up
close all;
clear variables;
clc;
format;


%% 
N = 1000;
alpha = 1; %3.983;
x = linspace(0,10,N);
Z = 1;
n = 1;
l = 0;
m = 0;

norm = sqrt(4)*(alpha)^(3/2); %((2*Z)/n)^(1/4) * sqrt(factorial(n-l-1)/(2*n*factorial(n+l)));

s1 = @(x) norm .* exp(-alpha.*x);
y = s1(x);


nor = @(a) sqrt(4*pi)*(2*a/pi).^(3./4);
G   = @(c,a,x) nor(a) * c * exp(-a.*x.^2); 

g4_2 = @(x) G(0.4301280,11.5356690,x)       + G(0.6789140,2.0533430,x);
g4_3 = @(x) G(0.15432897,30.1678710,x)      + G(0.53532814,5.4951153,x)       + G(0.44463454,1.4871927,x);
g4_6 = @(x) G(0.00916359628,312.87049370,x) + G(0.04936149294,57.364462530,x) + G(0.16853830490,16.048509400,x) + G(0.37056279970,5.5130961190,x) + G(0.41649152980,2.1408965530,x) + G(0.13033408410,0.8817394283,x);

g1_2 = @(x) G(0.430128498,1.309756377,x)    + G(0.678913531,0.233135974,x);
g1_3 = @(x) G(0.15432897,3.42525091,x)      + G(0.53532814,0.62391373,x)      + G(0.44463454,0.16885540,x);
g1_6 = @(x) G(0.00916359628,35.52322122,x)  + G(0.04936149294,6.513143725,x)  + G(0.16853830490,1.822142904,x) + G(0.37056279970,0.625955266,x) + G(0.41649152980,0.243076747,x) + G(0.13033408410,0.100112428,x);

dx = mean(diff(x));
Is1   = sum(dx * x.^2.*s1  (x).^2)
Ig4_2 = sum(dx * x.^2.*g4_2(x).^2)
Ig4_3 = sum(dx * x.^2.*g4_3(x).^2)
Ig4_6 = sum(dx * x.^2.*g4_6(x).^2)

Ig1_2 = sum(dx * x.^2.*g1_2(x).^2)
Ig1_3 = sum(dx * x.^2.*g1_3(x).^2)
Ig1_6 = sum(dx * x.^2.*g1_6(x).^2)

plot(x,x.^2.*s1(x),'k--', 'DisplayName', '1s')
grid('on')
hold('on')

%plot(x,x.^2.*g4_2(x)./Ig4_2, 'DisplayName', 'STO-2G')
%plot(x,x.^2.*g4_3(x)./Ig4_3, 'DisplayName', 'STO-3G')
%plot(x,x.^2.*g4_6(x)./Ig4_6, 'DisplayName', 'STO-6G')
plot(x,x.^2.*g1_2(x)./Ig1_2, 'DisplayName', 'STO-2G')
plot(x,x.^2.*g1_3(x)./Ig1_3, 'DisplayName', 'STO-3G')
plot(x,x.^2.*g1_3(x)./Ig1_6, 'DisplayName', 'STO-6G')

%h = legend('Hydrogen','STO-3G','STO-6G');
%h = legend('Hydrogen','STO-2G');
%set(h, 'FontSize', 16,'interpreter','latex');
h = legend('show');
set(h,'FontSize', 18, 'interpreter', 'latex');






