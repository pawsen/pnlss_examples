% Model input/output data of a Hammerstein system 
% from file fLMnlssWeighted.m

srcdir = '~/src/matlab';
addpath(genpath([srcdir '/PNLSS_v1_0']))

N = 2e3; % Number of samples
rng(10)
u = randn(N,1); % Input signal
f_NL = @(x) x + 0.2*x.^2 + 0.1*x.^3; % Nonlinear function
[b,a] = cheby1(2,5,2*0.3); % Filter coefficients
x = f_NL(u); % Intermediate signal
y = filter(b,a,x); % Output signal
scale = u\x; % Scale factor
sys = ss(tf(scale*b,a,[])); % Initial linear model = scale factor times underlying dynamics
nx = [2 3]; % 
ny = [2 3]; % Quadratic and cubic terms in output equation
T1 = 0; % No periodic signal transient handling
T2 = 200; % Number of transient samples to discard
model = fCreateNLSSmodel(sys.a,sys.b,sys.c,sys.d,nx,ny,T1,T2); % Initial linear model
model.xactive = fSelectActive('inputsonly',2,1,2,nx); % A Hammerstein system only has nonlinear terms in the input
model.yactive = fSelectActive('inputsonly',2,1,1,nx); % A Hammerstein system only has nonlinear terms in the input
MaxCount = 50; % Maximum number of iterations
W = []; % No weighting
[modelOpt,yOpt] = fLMnlssWeighted(u,y,model,MaxCount,W); % Optimized model and modeled output

rmse = @(y, yhat) sqrt(mean((y - yhat).^2));
fprintf('RMS error: %e\n', rmse(y, yOpt))

t = 0:N-1;
figure
plot(t,y,'b')
hold on
plot(t,yOpt,'r--')
plot(t,(y-yOpt),'k','LineWidth',3)
xlabel('Time')
ylabel('Output')
legend('True','Modeled','error')