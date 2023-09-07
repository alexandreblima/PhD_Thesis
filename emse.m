%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
%  Alexandre B. de Lima                                 %
%                                                       %
%   04/11/2007                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function emse() 
%
% Given a time series FARIMA(p,d,q), an estimate of the long memory
% parameter d (= H-1/2), and an adjusted AR(p2), this function returns: 
%   1) "k-steps-ahead-forecasts"  based on the last 100 observations
%      (k=1,...,k_max=100) given by the AR(p2) and the TARFIMA(L) models.
%   2) an Excel Table with the Empirical Mean-Square Errors (EMSE) for
%      several values of k_max (k_max= 7,10,20,50,100)
%
% Note that TARFIMA(L) is a finite-dimensional (i.e.,truncated) 
% state-space representation of an FARIMA(p,d,0) model. For more details, 
% please refer to the paper  
% "State-Space Modeling of Long-Range Dependent Teletraffic", ITC 2007,
% LNCS 4516, pp.260-271, 2007 (20th Int'l Teletraffic Congress, Ottawa, 
% Canada).
%
% inputs: 
%          timeseries = name of the text file which contains the time 
%                       series under analysis (ASCII column vector)
%          d          = estimated value of "d" for "timeseries" (I use 
%                       the estimate given by S+FinMetrics) 
%          name_ar_farima = name of a text file with the estimated AR coef.
%                       of FARIMA (if applicable)
%          name_ma_farima = name of a text file with the estimated MA coef.
%                       of FARIMA (if applicable)
%          L          = order of the finite-dimensional representation
%          namear     = name of a text file which contains the parameters 
%                       of the AR fitted model 
%
% outputs:
%          Table with EMSEs, Figures with forecasts 
%
% Copyright(c) 2007 Universidade de São Paulo, Laboratório de Comunicações
% e Sinais
% 
% Permission is granted for use and non-profit distribution providing that 
% this notice be clearly maintained. The right to distribute any portion 
% for profit or as part of any commercial product is specifically reserved
% for the author.

%function emse() 

clear
close all

% note that : a) vectors with coefficients are column vectors, and
%             b) vetors with observations are row vectors



% =========================================================================
% INPUT DATA

%  Prompt for time series
timeseries = input('file with the time series: ','s');

%  Prompt for estimated FARIMA model
d = input('parameter "d": ');

yes = 'y';
yes2 = 'Y';
j1 = input('AR coefficients? Y/N [Y]: ','s');
if ( strcmp(yes,j1) | strcmp(yes2,j1) ) 
    name_ar_farima = input('file with the estimated AR coef. of FARIMA: ','s');
    ar_farima = load(name_ar_farima);
    ar_farima = -1*ar_farima; % remember that the coefficients given by S-PLUS are the 
                          % negative of the coefficients used in MATLAB (convention)
                          % See S-PLUS Guide to Statistics, Vol. 2 
    ar_farima = [1; ar_farima]; % MATLAB convention
else
    ar_farima = [1];
end;

j2=input('MA coefficients? Y/N [Y]: ','s');
if ( strcmp(yes,j2) | strcmp(yes2,j2) )
    name_ma_farima = input('file with the estimated MA coef. of FARIMA: ','s');
    ma_farima = load(name_ma_farima);
    ma_farima = -1*ma_farima; % remember that the coefficients given by S-PLUS are the 
                          % negative of the coefficients used in MATLAB (convention)
    ma_farima = [1; ma_farima]; % MATLAB convention
else 
    ma_farima = [1];
end;

L = input('order of the finite-dimensional representation: ');
namear = input('file with the coef. of the estimated AR model: ','s');

fullseries =load(timeseries)';
N = length(fullseries);
%mu = mean(fullseries);

% =========================================================================
% Prediction from T-FARIMA(L) model

for j = 0:100 
    D(j+1) = -gamma(j-d)/(gamma(j+1)*gamma(-d));
end

D = -1*D; % D(B) = (1-B)^d = 1 -d_1B -d_2B^2 - ...
% "D" is a row vector

phi_temp = conv(ar_farima,D');

% Truncated AR(L) polynomial (1-phi_1B-phi_2B^2 - ... - phi_LB^L)
phi_B = ldiv(phi_temp,ma_farima,(L+1)); 
phi_B2 = -1*phi_B(2:(L+1)); %coef. of AR(L) (S+ convention) 
                            %column vector
[m,n]=size(phi_B2);
 
if (m<n)
    phi_B2 = phi_B2';
end

k = [7 10 20 50 100];
% k denotes the number of steps to predict ahead
% k_max=100

% ma100 = ldiv(ma_farima,phi_temp,101); % ldiv(b,a,N) is a function for 
                                        % the power series expansion of 
                                        % polynomial phi_B
% Now we must invert and truncate the AR(L) representation
% Yields the same result of: ma100 = ldiv(ma_farima,phi_temp,101)
ma100 = ldiv(1,phi_B,101);
ma100 = ma100(2:length(ma100));
% Vector with the coef. of the MA(100) representation of T-FARIMA(L) 
% [psi1 psi2 ... psi100], calculated by S-PLUS
% do not change signs! (as stated in section "DETAILS" of the "arma2ma()"
%                       function help )

%prev = zeros(length(fullseries)); % initialize row vector with predictions 

ma100quad = ma100.^2; % coefficiente to the square
var_erro = zeros(1,k(1)); % initialize row vector with the variances 
                          % of the k-step ahead prediction errors

%t = 2700;
%t=866;
% origin of prediction
                          
for j=1:length(k)

    t=length(fullseries)-k(j); 
    % origin of prediction
    
    
    % estimate mean of time series
    mu = mean(fullseries(t:-1:1));  
    
    % Initialize row vector Z (with L previous observations)
    Z  = fullseries(t:-1:t-(length(phi_B2)-1)) - mu;
    % Z = [Z_358 Z_357 ... Z_309]

    % column vector I of innovations
    I = filter(phi_B,1,(fullseries'-mu));
    pot_inov = var(I);

    for h=1:k(j)    
        prev(h) = Z*phi_B2;
        Z = [prev(h) Z(1:(length(phi_B2)-1))];
    end    

    prev = prev + mu; % predictions
    x = fullseries(t+1:t+k(j)); % future observations
    EQMPt(j) = sum((prev-x).^2)/k(j);
    erro = x-prev;
    
    if j == length(k)
        for h=1:k(j)    
            if h==1
                var_erro(h) = 1*pot_inov; % vide Eq.(9.11) Morettin(2004)
            else
                var_erro(h) = (ma100quad(h-1)*pot_inov) + var_erro(h-1);
            end 
        end   
    end   
    
end

upper_bound = prev + 1.96*sqrt(var_erro);
lower_bound = prev - 1.96*sqrt(var_erro);

erroabs = abs(erro);

% =========================================================================

% Prediction from AR(p) model 

ar = load(namear);
p = length(ar);
ar = -1*ar; % remember that the coefficients given by S-PLUS are the 
            % negative of the coefficients used in MATLAB (convention)
            
ar = [1; ar]; % MATLAB convention
ARorder = strcat('AR(',num2str(p),')')
B = [1];
[z_ar,p_ar,G] = tf2zpk(B,ar);
zplane(z_ar,p_ar,'k')
title(['Diagrama de Pólos e Zeros do modelo AR(',int2str(p),')'])
figure

ar2 = -1*ar(2:length(ar)); % S-PLUS convention

ar_ma100 = ldiv(1,ar,101);
ar_ma100 = ar_ma100(2:length(ar_ma100))';
% Vector with the coef. of the MA(100) representation of AR(p) 
% [psi1 psi2 ... psi100], calculated by S-PLUS
% do not change signs! (as stated in section "DETAILS" of the "arma2ma()"
%                       function help )

%prev2 = zeros(length(fullseries)); % initialize row vector with predictions 
ar_ma100quad = ar_ma100.^2; % coef. to the square
var_erro2 = zeros(1,k(1));  % initialize row vector with the variances 
                            % of the k-step ahead prediction errors


for j=1:length(k)
    
    t=length(fullseries)-k(j); % origin of prediction
    
    % estimate mean of time series
    mu = mean(fullseries(t:-1:1));  
    
    % Initialize row vector Z2 (p observations)
    Z2 = fullseries(t:-1:t-(length(ar2)-1))- mu;

    % column vector I2 of innovations
    I2 = filter(ar,1,(fullseries'-mu));
    pot_inov2 = var(I2);
    
    for h=1:1:k(j)    
        prev2(h) = Z2*ar2;
        Z2 = [prev2(h) Z2(1:(length(ar2)-1))];
    end

    prev2 = prev2 + mu;
    x2 = fullseries(t+1:t+k(j));
    EQMPt2(j) = sum((prev2-x2).^2)/k(j);  
    erro2 = x2-prev2;

    if j == length(k)
        for h=1:k(j)    
            if h==1
                var_erro2(h) = 1*pot_inov2; % vide Eq.(9.11) Morettin(2004)
            else
                var_erro2(h) = (ar_ma100quad(h-1)*pot_inov2) + var_erro2(h-1);
            end 
        end   
    end   
    
end

upper_bound2 = prev2 + 1.96*sqrt(var_erro2);
lower_bound2 = prev2 - 1.96*sqrt(var_erro2);

erro2abs = abs(erro2);

delta_error_abs = erro2abs - erroabs ;
% =========================================================================
%ARorder = strcat('AR(',num2str(p),')')
TFARIMAorder = strcat('TARFIMA(',num2str(L),')')

EQM = [EQMPt2',EQMPt'];

emse1 = {'k',ARorder,TFARIMAorder; k(1),EQM(1,1),EQM(1,2); k(2),EQM(2,1),EQM(2,2); k(3),EQM(3,1),EQM(3,2);
        k(4),EQM(4,1),EQM(4,2); k(5),EQM(5,1),EQM(5,2)}

xlswrite('emse.xls', emse1, 'teste');

%=========================================================================
% Draw graphics

% Estimation of Power spectral density (PSD) using Welch's method
% Advantages: simultaneous reduction of a) bias due to leakage (via data 
% tapering) and b) variability (by averaging the periodograms of 
% Nb blocks of size Ns) of the periodogram.

[Pxx,w] = pwelch(fullseries);
plot(w/(2*pi),10*log10(Pxx))
%xlabel('normalized frequency (cycle/sample)')
xlabel('f (ciclo/amostra)')
%ylabel('power/frequency (dB/cycle/sample)')
ylabel('dB/ciclo/amostra')
%title('PSD via WOSA')
title('DEP estimada (método WOSA)')

figure
[Pxx,w] = pwelch(I);
plot(w/(2*pi),10*log10(Pxx))
xlabel('f (ciclo/amostra)')
ylabel('dB/ciclo/amostra')
title('DEP estimada (método WOSA) dos resíduos do modelo TARFIMA')

t=length(fullseries)-max(k); % origin of prediction using k_max=100
                             % t=1106 for the Nile River minima series
                             % used by Beran in the first chapter of
                             % his book

                             
%n = (t-max(k)+1):1:N;
n = (t-max(k)+1):1:(t+max(k));

figure
plot(n,fullseries(n),'k.-')
hold on
%pause
%plot((t+1:N),prev,'.:g')
plot((t+1:t+max(k)),prev,'.-g')
hold on
%pause
%plot((t+1:N),upper_bound,'x:g')
plot((t+1:t+max(k)),upper_bound,'x:g')
hold on
%pause
%plot((t+1:N),lower_bound,'*:g')
plot((t+1:t+max(k)),lower_bound,'*:g')
hold on
%pause
%plot((t+1:N),prev2,'.:r')
plot((t+1:t+max(k)),prev2,'.:r')
hold on
%pause
%plot((t+1:N),upper_bound2,'x:r')
plot((t+1:t+max(k)),upper_bound2,'x:r')
hold on
%pause
%plot((t+1:N),lower_bound2,'*:r')
plot((t+1:t+max(k)),lower_bound2,'*:r')
legend('observação',['previsão TARFIMA(',int2str(L),')'], ['limite superior TARFIMA(',int2str(L),')'],['limite inferior TARFIMA(',int2str(L),')'],['previsão AR(',int2str(p),')'],['limite superior AR(',int2str(p),')'],['limite inferior AR(',int2str(p),')'])
%legend('observations',['TARFIMA(',int2str(L),') forecasts'], ['upper limit TARFIMA(',int2str(L),')'],['lower limit TARFIMA(',int2str(L),')'],['AR(',int2str(p),') forecasts'],['upper limit AR(',int2str(p),')'],['lower limit AR(',int2str(p),')'])


figure
plot((t+1:t+max(k)),delta_error_abs,'k.-')
xlabel('h')
%ylabel('difference of absolute values of prediction errors')
ylabel('diferença entre os erros absolutos de previsão')

%=========================================================================
% Comparison of TARFIMAs with L=10, and 100

L2 = 10;
phi_B3 = ldiv(phi_temp,ma_farima,(L2+1)); 
phi_B4 = -1*phi_B3(2:(L2+1)); 
[m,n]=size(phi_B4);
 
if (m<n)
    phi_B4 = phi_B4';
end

ma100L2 = ldiv(1,phi_B3,101);
ma100L2 = ma100L2(2:length(ma100L2));

ma100quadL2 = ma100L2.^2; 
var_erro3 = zeros(1,max(k)); 

mu = mean(fullseries(t:-1:1));  
Z3  = fullseries(t:-1:t-(length(phi_B4)-1)) - mu;
I3 = filter(phi_B3,1,(fullseries'-mu));
pot_inov3 = var(I3);

for h=1:max(k)    
    prev3(h) = Z3*phi_B4;
    Z3 = [prev3(h) Z3(1:(length(phi_B4)-1))];
end    

prev3 = prev3 + mu; % predictions
x = fullseries(t+1:t+max(k)); % future observations
    
erro3 = x-prev3;
    
for h=1:max(k)    
    if h==1
       var_erro3(h) = 1*pot_inov3; % vide Eq.(9.11) Morettin(2004)
    else
       var_erro3(h) = (ma100quadL2(h-1)*pot_inov2) + var_erro3(h-1);
    end 
end   
                          
upper_bound3 = prev3 + 1.96*sqrt(var_erro3);
lower_bound3 = prev3 - 1.96*sqrt(var_erro3);         

n = (t-max(k)+1):1:(t+max(k));

figure
plot(n,fullseries(n),'k.-')
hold on
%pause
%plot((t+1:N),prev,'.:g')
plot((t+1:t+max(k)),prev,'.:g')
hold on
%pause
%plot((t+1:N),upper_bound,'x:g')
plot((t+1:t+max(k)),upper_bound,'x:g')
hold on
%pause
%plot((t+1:N),lower_bound,'*:g')
plot((t+1:t+max(k)),lower_bound,'*:g')
hold on
%pause
%plot((t+1:N),prev2,'.:r')
plot((t+1:t+max(k)),prev3,'.:b')
hold on
%pause
%plot((t+1:N),upper_bound2,'x:r')
plot((t+1:t+max(k)),upper_bound3,'x:b')
hold on
%pause
%plot((t+1:N),lower_bound2,'*:r')
plot((t+1:t+max(k)),lower_bound3,'*:b')
%legend('observations',['TARFIMA(',int2str(L),') forecasts'], ['upper limit TARFIMA(',int2str(L),')'],['lower limit TARFIMA(',int2str(L),')'],'TARFIMA(10)','upper limit TARFIMA(10)','lower limit TARFIMA(10)')
legend('observação',['previsão TARFIMA(',int2str(L),')'], ['limite superior TARFIMA(',int2str(L),')'],['limite inferior TARFIMA(',int2str(L),')'],'previsão TARFIMA(10)','limite superior TARFIMA(10)','limite inferior TARFIMA(10)')


