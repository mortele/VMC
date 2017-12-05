%% Clean up
close all;
%clear variables;
clc;
format;

%% Parameters
n = 2;
N = 2500;
Z = 4;
a = 3.983; %1.843;%10.22;%3.983;
alpha = a;
x1 = logspace(0,log10(4.5+1),N)-1;
x1 = x1';
x2 = logspace(0,log10(6+1),N)-1;
x2 = x2';
x3 = logspace(0,log10(2+1),N)-1;
x3 = x3';
A  = logical([ones(1,n) zeros(1,n)]);
C  = logical([zeros(1,n) ones(1,n)]);

% Limits on the coefficients.
lowerC = -inf;
upperC =  inf;
lowerA =    0;
upperA =  inf;
lower = @(n) [ones(1,n)*lowerA ones(1,n)*lowerC];
upper = @(n) [ones(1,n)*upperA ones(1,n)*upperC];


%% Wave function
% Node-less Slater type orbitals
norm1S = sqrt(a^3/pi);
norm2S = (1/4.)*sqrt(a^5/(6*pi));
norm2P = (1/8.)*sqrt(a^7/(15*pi));
s1 = @(x) norm1S .* exp(-a.*x);
s2 = @(x) norm2S .* x .* exp(-a * 0.5 .* x);
p2 = @(x) norm2P .* x .* exp(-a * 0.5 .* x);

% Hydrogen orbitals
% norm1S = sqrt(a^3/pi);
% norm2S = sqrt(a^3/(8*pi));
% norm2P = sqrt(a^5/(32*pi));
% s1 = @(x) norm1S .* exp(-a.*x);
% s2 = @(x) norm2S .* (1-a*0.5*x) .* exp(-a * 0.5 .* x);
% p2 = @(x) norm2P .* (1-a*0.5*x) .* exp(-a * 0.5 .* x);

y1 = s1(x1);
y2 = s2(x2);
y3 = p2(x3);

% Prepare figures for plotting.
figure(1);
plot(x1,y1,'k--','DisplayName','1s');
hold on;
if Z > 3
    figure(2);
    plot(x2,y2,'k--','DisplayName','2s');
    hold on;
end
if Z > 5 
    figure(3);
    plot(x3,y3,'k--','DisplayName','2p');
    hold on;
end

%% Primitive function
g = @(c,a,x) c*exp(-a*x.^2);


%% STO-nG
% 1s orbital
functionalForm = {'c1*exp(-a1*x.^2)',
    'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)',
    'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)+c3*exp(-a3*x.^2)',
    'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)+c3*exp(-a3*x.^2)+c4*exp(-a4*x.^2)',
    'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)+c3*exp(-a3*x.^2)+c4*exp(-a4*x.^2)+c5*exp(-a5*x.^2)',
    'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)+c3*exp(-a3*x.^2)+c4*exp(-a4*x.^2)+c5*exp(-a5*x.^2)+c6*exp(-a6*x.^2)',
    'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)+c3*exp(-a3*x.^2)+c4*exp(-a4*x.^2)+c5*exp(-a5*x.^2)+c6*exp(-a6*x.^2)+c7*exp(-a7*x.^2)',
    'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)+c3*exp(-a3*x.^2)+c4*exp(-a4*x.^2)+c5*exp(-a5*x.^2)+c6*exp(-a6*x.^2)+c7*exp(-a7*x.^2)+c8*exp(-a8*x.^2)',
    'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)+c3*exp(-a3*x.^2)+c4*exp(-a4*x.^2)+c5*exp(-a5*x.^2)+c6*exp(-a6*x.^2)+c7*exp(-a7*x.^2)+c8*exp(-a8*x.^2)+c9*exp(-a9*x.^2)',
    'c1*exp(-a1*x.^2)+c2*exp(-a2*x.^2)+c3*exp(-a3*x.^2)+c4*exp(-a4*x.^2)+c5*exp(-a5*x.^2)+c6*exp(-a6*x.^2)+c7*exp(-a7*x.^2)+c8*exp(-a8*x.^2)+c9*exp(-a9*x.^2)+c10*exp(-a10*x.^2)'};
%%
[f,gof] = fit(x1,y1,functionalForm{n},...
    'StartPoint',2*randn([2*n 1])-1,...
    'Lower',lower(n),...
    'Upper',upper(n),...
    'Algorithm','Trust-Region',...
    'Normalize','off',...
    'MaxIter',1e5,...
    'MaxFunEvals',1e5,...
    'DiffMinChange',1e-8,...
    'DiffMaxChange',0.1,...
    'Robust','LAR',...
    'TolFun',1e-6,...
    'TolX',1e-6)
clear STO1;
ca1 = coeffvalues(f);
a1  = ca1(A);
c1  = ca1(C);
STO1= g(c1(1),a1(1),x1);
for i=2:n
    STO1 = STO1 + g(c1(i),a1(i),x1);
end
figure(1);
plot(x1,STO1,'DisplayName',strcat('STO-',num2str(n),'G'));

c = c1;
a = a1;

%%
if Z > 3
    % 2s orbital
    f = fit(x2,y2,functionalForm{n},...
        'StartPoint',2*randn([2*n 1])-1,...
        'Lower',lower(n),...
        'Upper',upper(n),...
        'Algorithm','Trust-Region',...
        'Normalize','off',...
        'MaxIter',1e5,...
        'MaxFunEvals',1e5,...
        'DiffMinChange',1e-8,...
        'DiffMaxChange',0.1,...
        'Robust','LAR',...
        'TolFun',1e-6,...
        'TolX',1e-6)
    ca2 = coeffvalues(f);
    a2  = ca2(A);
    c2  = ca2(C);
    STO2= g(c2(1),a2(1),x2);
    for i=2:n
        STO2 = STO2 + g(c2(i),a2(i),x2);
    end
    
    figure(2);
    plot(x2,STO2,'DisplayName',strcat('STO-',num2str(n),'G'));
    
    
    c = [c1;c2];
    a = [a1;a2];
end


%%
if Z > 5
    % 2p orbital
    ca3 = coeffvalues(f);
    a3  = ca3(A);
    c3  = ca3(C);
    
    % No need to fit another function, we can simply take the old
    % coefficients for 2S and rescale them according to the 2P
    % normalization constant.
    c3  = c3 .* (norm2P/norm2S);
    STO3= g(c3(1),a3(1),x3);
    for i=2:n
        STO3 = STO3 + g(c3(i),a3(i),x3);
    end
    
    figure(3);
    plot(x3,STO3,'DisplayName',strcat('STO-',num2str(n),'G'));
    
    
    c = [c1;c2;c3];
    a = [a1;a2;a3];
end


%%
%writeBasisToFile(Z,n,c,a);


%% Plot esthetics
ee = 1;
if Z > 3
    ee = 2;
end
if Z > 5
    ee = 3;
end
for i=1:ee
    figure(i);
    h = legend('show');
    set(h,'FontSize',18,'interpreter','latex');
    xlabel('$r [a_0]$','FontSize',18,'interpreter','latex');
    ylabel('Amplitude','FontSize',18,'interpreter','latex');
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 18)
    xt = get(gca, 'YTick');
    set(gca, 'FontSize', 18)
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    if i==1
        axis([0 3 0 1.05]);
    end
   
end


%%
% x1 = x1(100:end);
% x2 = x2(100:end);
% lap1s = norm1S .* exp(-alpha.*x1).*alpha.*(x1.*alpha-2)./x1;
% lap2s = norm2S .* exp(-alpha.*x2/2).*(8+x2*alpha.*(x2.*alpha-8))./(4.*x2);
% glap = @(c,a,x) 2.*exp(-a.*x.^2).*a.*(-3 + 2.*x.^2.*a);
% 
% glap1s = glap(c1(1),a1(1),x1);
% for i=2:n
%     glap1s = glap1s + glap(c1(i),a1(i),x1);
% end
% glap2s = glap(c2(1),a2(1),x2);
% for i=2:n
%     glap2s = glap2s + glap(c2(i),a2(i),x2);
% end
% 
% RMSE_lap_1s = sqrt(mean((glap1s-lap1s).^2))
% RMSE_lap_2s = sqrt(mean((glap2s-lap1s).^2))
% 
% figure;
% plot(x1,lap1s,'k--');
% hold on;
% plot(x1,glap1s,'r-');
% 
% figure;
% plot(x2,lap2s,'k--');
% hold on;
% plot(x2,glap2s,'r-');

%%
%close all;
figure(10);
hold on;
plot(x1,y1,'k--','DisplayName','1s');
%plot(x1,STO11,'DisplayName','STO-1G');
%plot(x1,STO22,'DisplayName','STO-2G');
%plot(x1,STO33,'DisplayName','STO-3G');
%plot(x1,STO44,'DisplayName','STO-4G');
plot(x1,STO55,'DisplayName','STO-5G');
plot(x1,STO66,'DisplayName','STO-6G');
plot(x1,STO77,'DisplayName','STO-7G');
%plot(x1,STO88,'DisplayName','STO-8G');

h = legend('show');
set(h,'FontSize',18,'interpreter','latex');
xlabel('$r [a_0]$','FontSize',18,'interpreter','latex');
ylabel('Amplitude','FontSize',18,'interpreter','latex');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 18)
xt = get(gca, 'YTick');
set(gca, 'FontSize', 18)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
axis([0 3 0 1.05]);

figure;
eSTO1 = mean(abs(y1-STO11));
eSTO2 = mean(abs(y1-STO22));
eSTO3 = mean(abs(y1-STO33));
eSTO4 = mean(abs(y1-STO44));
eSTO5 = mean(abs(y1-STO55));
eSTO6 = mean(abs(y1-STO66));
eSTO7 = mean(abs(y1-STO77));


eee = [ eSTO1;
        eSTO2;
        eSTO3;
        eSTO4;
        eSTO5;
        eSTO6;
        eSTO7];
ppp = [1 2 3 4 5 6 7];

semilogy(ppp,eee,'ko-');

%            sse: 5.3550e-05
%        rsquare: 1.0000
%            dfe: 2488                     6
%     adjrsquare: 1.0000
%           rmse: 1.4671e-04

%            sse: 1.7059e-05
%        rsquare: 1.0000
%            dfe: 2486                     7
%     adjrsquare: 1.0000 
%           rmse: 8.2838e-05


%h = legend('show');
%set(h,'FontSize',18,'interpreter','latex');
xlabel('Primitives','FontSize',18,'interpreter','latex');
ylabel('mean absolute error w.r.t. 1s','FontSize',18,'interpreter','latex');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 18)
xt = get(gca, 'YTick');
set(gca, 'FontSize', 18)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];


