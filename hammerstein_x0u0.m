% Model input/output data of a Hammerstein system with non-zero initial conditions
% from file fLMnlssWeighted_x0u0.m

srcdir = '~/src/matlab';
addpath(genpath([srcdir '/PNLSS_v1_0']))

N = 200; % Number of samples
NTrans = 100; % Number of samples after zero initial conditions
rng(10)
u = randn(NTrans+N,1); % Input signal
f_NL = @(x) x + 0.2*x.^2 + 0.1*x.^3; % Nonlinear function
[b,a] = cheby1(2,5,2*0.3); % Filter coefficients
x = f_NL(u); % Intermediate signal
y = filter(b,a,x); % Output signal
u(1:NTrans) = []; % Remove first NTrans samples to obtain non-zero intial conditions
x(1:NTrans) = []; % Remove first NTrans samples to obtain non-zero intial conditions
y(1:NTrans) = []; % Remove first NTrans samples to obtain non-zero intial conditions
scale = u\x; % Scale factor
sys = ss(tf(scale*b,a,[])); % Initial linear model = scale factor times underlying dynamics
nx = [2 3]; % Quadratic and cubic terms in state equation
ny = [2 3]; % Quadratic and cubic terms in output equation
T1 = 0; % No periodic signal transient handling
T2 = []; % No transient samples to discard
model = fCreateNLSSmodel(sys.a,sys.b,sys.c,sys.d,nx,ny,T1,T2); % Initial linear model
model.xactive = fSelectActive('inputsonly',2,1,2,nx); % A Hammerstein system only has nonlinear terms in the input
model.yactive = fSelectActive('inputsonly',2,1,1,nx); % A Hammerstein system only has nonlinear terms in the input
MaxCount = 50; % Maximum number of iterations
W = []; % No weighting
[modelOpt,yOpt] = fLMnlssWeighted(u,y,model,MaxCount,W); % Optimized model and modeled output (without estimating initial conditions)
model_x0u0 = model; % Estimate initial conditions
model_x0u0.x0active = (1:model_x0u0.n).'; % Estimate initial conditions
model_x0u0.u0active = (1:model_x0u0.m).'; % Estimate initial conditions
[modelOpt_x0u0,yOpt_x0u0] = fLMnlssWeighted_x0u0(u,y,model_x0u0,MaxCount,W); % Optimized model and modeled output (initial conditions estimated);

rmse = @(y, yhat) sqrt(mean((y - yhat).^2));
fprintf('RMS error: %e\nRMS error(initial conditions estimated): %e\n',...
    rmse(y,yOpt), rmse(y,yOpt_x0u0))

t = 0:N-1;
figure
plot(t,y,'b')
hold on
plot(t,y-yOpt,'-r','LineWidth',3)
plot(t,y-yOpt_x0u0,'.g','LineWidth',3)
xlabel('Time')
ylabel('Output / output error')
legend('True','Error PNLSS','Error PNLSS (initial conditions estimated)')
