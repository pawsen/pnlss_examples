% -*- coding: utf-8 -*-
%
% Example of nonparametric system characterization
%
% This script estimates
% - Nonlinear distortion compared to noise level
% - Type of nonlinearity (odd or even)
%
% To do that, we excite the system using a multisine signal. To estimate the
% noise, we excite with P periods to obtain steady state and average over the
% steady state periods. To estimate nonlinear distortion, we excite with R
% realizations(distinct experiments) and average over the realizations.
% We average over the FRF, G = Y / U (also called BLA; the Best Linear Approx)
%
% See:
% Schoukens, J., Vaes, M., & Pintelon, R., Linear System Identification in a Nonlinear Setting:
% Nonparametric Analysis of the Nonlinear Distortions and Their Impact on the Best Linear Approximation, , 36(3), 38–69 ().
% http://dx.doi.org/10.1109/MCS.2016.2535918
%
% See also http://homepages.vub.ac.be/~ktiels/pnlss.html for a more
% robust implementation of the BLA calculation. For mechanical problems
% the simple way of just taking averages directly is good enough.

clc
clear variables

%% Settings for signal
fs = 1500;  % Hz
N = 4096;   % number of points per period
P = 3;      % number of periods
R = 3;     % number of realizations
fMin = 0;   % lowest excited frequency
fMax = 500; % highest excited frequency

uStd = 2;   % standard deviation of the generated signal
nStd = 0.01; % standard deviation of the noise added on the output

time = 0:1/fs:N*P/fs - 1/fs;
freq = 0:fs/N:fs-fs/N;

%% generate data
options.N = N;
options.P = P;
options.M = R;
options.fMin = fMin;   % Hz - first frequency above zere will be selected
options.fMax = fMax;   % Hz
options.fs = fs;       % Hz
options.type = 'full'; %'full', 'odd', 'oddrandom'
options.type = 'oddrandom';
options.std = uStd;
[u, lines] = fMultiSinGen(options);

%% Define system
% generate LTI block
[b,a] = cheby2(3,30,0.2);

% static nonlinearity at the input
x = zeros(size(u));
x = tanh(u); % smooth saturation NL
%x = u + 0.005*u.^2 - 0.01*u.^3; % a polynomial NL
%x = u  - 0.01*u.^3; 

% LTI filtering
y = filter(b,a,x);
y = y + nStd*randn(size(y));  % noise

%% estimate BLA

% reshape and remove 1 period because of transients
data.u = reshape(u(N+1:end,:),[N P-1 R]);
data.y = reshape(y(N+1:end,:),[N P-1 R]);

[BlaTemp] = fBlaEst(data,[]);
BlaEst = BlaTemp;

%% plot
% Thes two plots shows the same. This one uses G = Y/U, the next Y = G*U
figure; hold on;
plot(freq(lines),db(BlaEst.FRF(lines)))
plot(freq(lines),db(BlaEst.tVar(lines)*R,'power'),'s')
plot(freq(lines),db(BlaEst.nVar(lines)*R,'power'),'k-')
xlabel('frequency (Hz)')
ylabel('magnitude (dB)')
title(['Estimated BLA: uStd = ' num2str(uStd)])
legend('BLA FRF','Total Distortion','Noise Distortion')
% total and noise distortion averaged over P periods and R realizations
% total distortion level includes nonlinear and noise distortion

U = fft(data.u(:,end,1));
Y = fft(data.y(:,end,1));
figure; hold on;
plot(freq(lines),db(Y(lines)))
plot(freq(lines),db(BlaEst.FRF(lines).*U(lines)))
plot(freq(lines),db(R*BlaEst.tVar(lines).*abs(U(lines)).^2,'power'),'s')
plot(freq(lines),db(R*BlaEst.nVar(lines).*abs(U(lines)).^2,'power'),'k-')
xlabel('frequency (Hz)')
ylabel('magnitude (dB)')
title(['Estimated BLA Output: uStd = ' num2str(uStd)])
legend('Y','Y_{BLA}','Total Distortion Output','Noise Distortion Output')
% total and noise distortion averaged over P periods and 1 realizations
% total distortion level includes nonlinear and noise distortion

%%
% Distinguis odd and even NL
% Calculate 'not-excited' lines
non_lines = 1:N/2;
non_lines(lines) = [];
non_odd = non_lines(logical(mod(non_lines-1,2)));
non_even = non_lines(~mod(non_lines-1,2));
Y = fft(data.y,[],1);

% Plot the importance of the odd and even nonlinearity
figure; hold on;
plot(freq(lines),db(Y(lines)))
plot(freq(non_odd),db(Y(non_odd)))
plot(freq(non_even),db(Y(non_even)))
plot(freq(lines),db(R*BlaEst.nVar(lines).*abs(U(lines)).^2,'power'),'k')
legend('Y', 'non_odd', 'non_even', 'noise')


function [Bla] = fBlaEst(data)
% 
% estimates the BLA of a nonlinear system in the open loop case
%
% INPUT:
% data.u: NxPxM containing the input (N: # points per period, P: # periods, M: # realizations)
% data.y: NxPxM containing the transient free output
%
% OUTPUT:
% Bla.FRF: nLinesx1 containing the estimated FRF
% BLA.nVar: nLinesx1 containing the estimated noise variance
% BLA.tVar: nLinesx1 containing the estimated total variance (noise + variance)
% all variances are scaled with respect to one period and one realization

u = data.u;
y = data.y;

[N,P,M] = size(u);

U = fft(u,[],1);
Y = fft(y,[],1);

G = Y./U; % calculate BLA for every realization adn period
Gm = squeeze(mean(G,2)); % average over the periods
Gbla = squeeze(mean(Gm,2)); % average over the realizations

tVar = var(Gm,[],2); % variance over the realizatioin: total variance
nVar = mean(squeeze(var(G,[],2)),2); % average variance over the periods: noise variance

Bla.FRF = Gbla;
Bla.nVar = nVar/M/P;
Bla.tVar = tVar/M;
end