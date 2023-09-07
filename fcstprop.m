%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
%  Alexandre B. de Lima                                 %
%                                                       %
%   14/10/2007                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script gives measures of the k-step forecastability of a time series
% as defined in p.24 of "An Introduction to Long-Memory Time Series Models 
% and Fractional Differencing", Granger and Joyeux, 1980. 
% 
% Copyright(c) 2007 Universidade de São Paulo, Laboratório de Comunicações
% e Sinais
% 
% Permission is granted for use and non-profit distribution providing that 
% this notice be clearly maintained. The right to distribute any portion 
% for profit or as part of any commercial product is specifically reserved
% for the author.

clear
close all

% note that : a) vectors with coefficients are column vectors, and
%             b) vetors with observations are row vectors

% =========================================================================

% Forecasting potential of the long-memory models (predictions based on 
% an infinite past)

ar_farima = [1];
ma_farima = [1];
d = input('parameter "d": '); % parameter "d"

M = 200;  % to derive the first M coef. of the MA(inf) representation

for j = 0:M 
    D(j+1) = -gamma(j-d)/(gamma(j+1)*gamma(-d)); 
    % Eq.(16.30), pág. 475, "Análise de Séries Temporais", 
    % Morettin & Toloi,2004
end

D = -1*D; % D(B) = (1-B)^d = 1 -d_1B -d_2B^2 - ... (estimated truncated 
          % fractional difference filter)
          % "D" is a row vector 

phi_B = conv(ar_farima,D'); % AR polynomial (1-phi_1B-phi_2B^2 - ...) 
phi_B2 = -1*phi_B(2:length(phi_B)); %coef. of AR polynomial (S+ convention) 

ma100 = ldiv(ma_farima,phi_B,M); % ldiv(b,a,N) is a function for 
                                   % the power series expansion of 
                                   % polynomial phi_B
ma100 = ma100';
ma100quad = ma100.^2; % coefficients to the square

Sk = cumsum(ma100quad); % vector with the variances of the k-step forecast 
                        % errors (length = 100) for the ARFIMA model

% let x_t be an FARIMA(0,d,0) of the form 
%       (1-B)^d{x_t} = w_t
% where it is assumed that var{w_t}=1 for convenience.
% The variance of x_t is given by 
%       var{x_t} = gamma(-2d+1)/gamma^2(-d+1)
% (see, Eq.(3.2) in paper "Fractional Differencing", Hosking, 1981) 
Vd = gamma(-2*d+1)/(gamma(-d+1))^2 ;                       

k = [1 10 20 100]; % so that you can reproduce the results presented in 
                   % Tables 8.7 and 8.8 by Jan Beran in "Statistics for 
                   % Long-Memory Processes", Chapman & Hall, 1994

% quantity "Rk" measures the k-step ahead forecasting potential of the LRD
% model, i.e., it's a standardized measure of how well one can predict 
% x_{t+k} given the infinite past

for j=1:length(k)
    Rk(j) = 1 - Sk(k(j))/Vd; 
end

%==========================================================================
% Forecasting potential of the TARFIMA(L), whose predictions are based on
% L previous observations

% INPUT DATA

L = input('order of the finite-dimensional representation: ');
%d2 = input('estimate of "d": ');

for j = 0:L 
    D2(j+1) = -gamma(j-d)/(gamma(j+1)*gamma(-d));
end

D2 = -1*D2; 

phi2_B = conv(ar_farima,D2'); % AR polynomial (1-phi_1.B-phi_2.B^2 - ...) 
phi2_B2 = -1*phi2_B(2:length(phi2_B)); %coef. of AR polynomial (S+ convention) 

ma100_2 = ldiv(ma_farima,phi2_B,M); 
ma100_2 = ma100_2';
ma100quad_2 = ma100_2.^2;                                    
                                    
Sk2 = cumsum(ma100quad_2); 
% vector with the variances of the k-step forecast errors (length = 100)
% for the TARFIMA model  

% let x_t be approximated by an TARFIMA(L) model as 
%       (1 - phi_1.B - phi_2.B^2 - ... - phi_L.B^L)x_t = w_t
% where it is assumed that var{w_t}=1 for convenience.
% x_t has the MA(inf) representation
%       x_t = (1 - phi_1.B - phi_2.B^2 - ... - phi_L.B^L)^{-1}w_t  or
%       x_t = (1 + psi_1.B + psi_2.B^2 + ... )w_t
% and the variance of x_t is given by 
%       var{x_t} = 1 + sum_{j=1}^{inf}(psi_j^2)
% (see, Eq.(3.3.1) in p.91 of "Time Series: Theory and Methods", 2nd Ed., 
%  Brockwell and Davis, 1991)
% Below, Sk2(M) gives an estimate of var{x_t}
Vd2 =  Sk2(M);  

% quantity "R2k" measures the k-step ahead forecasting potential of the 
% truncated FARIMA model, i.e., it's a standardized measure of how well one 
% can predict x_{t+k} given "L" (=order of TARFIMA model) past observations
for j=1:length(k)
    R2k(j) = 1 - Sk2(k(j))/Vd2; 
end

% =========================================================================

% Forecasting potential for the AR(1) process with with rho(1)=1/9 
% and rho(1)=2/3. The lag-1 correlation is chosen such that it is the same
% as for a ARFIMA(0,d,0) process with H=0.6 and 0.9 respectively.
rho_1 = 1/9;
rho_2 = 2/3;

% let x_t be an AR(1) of the form 
%       (1-phi_1)^d{x_t} = w_t
% It is assumed that var{w_t}=1 for convenience.

Vh06_AR1 = 1/(1 -rho_1^2); % variance of AR(1) process with H=0.6;
Vh09_AR1 = 1/(1 -rho_2^2); % variance of AR(1) process with H=0.9;
                           % see Brockwell & Davis (1991), 
                           % Morettin & Toloi (2004)

for j=1:100                           
    Skh06_AR1(j) = Vh06_AR1*(1 - rho_1^(2*j));
    Skh09_AR1(j) = Vh09_AR1*(1 - rho_2^(2*j));
end

for j=1:length(k)
    R3k(j) = 1 - Skh06_AR1(k(j))/Vh06_AR1; 
    R4k(j) = 1 - Skh09_AR1(k(j))/Vh09_AR1;
end




