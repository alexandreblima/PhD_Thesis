%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
%  Alexandre B. de Lima                                 %
%                                                       %
%   28/10/2007                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

j=sqrt(-1);

d = input('parameter "d": ');
%d=0.3859;
H = d + 0.5;
da = 1; %potencia do ruido branco Gaussiano

yes = 'y';
yes2 = 'Y';
j1 = input('AR coefficients? Y/N [Y]: ','s');
if ( strcmp(yes,j1) | strcmp(yes2,j1) ) 
    name_ar_farima = input('file with the estimated AR coef. of FARIMA: ','s');
    ar_farima = load(name_ar_farima);
else
    ar_farima = [1];
end;

j2=input('MA coefficients? Y/N [Y]: ','s');
if ( strcmp(yes,j2) | strcmp(yes2,j2) )
    name_ma_farima = input('file with the estimated MA coef. of FARIMA: ','s');
    ma_farima = load(name_ma_farima);
else 
    ma_farima = [1];
end;

ar_farima = -1*ar_farima; % remember that the coefficients given by S-PLUS are the 
                          % negative of the coefficients used in MATLAB (convention)
                          % See S-PLUS Guide to Statistics, Vol. 2 
ar_farima = [1; ar_farima];
ma_farima = -1*ma_farima;                          
ma_farima = [1; ma_farima];

N=2048;
Deltaf = (0.5/N);

i=0;

% Calculo da Densidade Espectral de Potencia do modelo FARIMA(p,d,0) estimado
for f=(0.5/N):(0.5/N):0.5
    A = [1 -exp(-j*2*pi*f)];
    B = (abs(sum(A))^(-2*d)*(da^2));
    if (length(ma_farima)>1)
       A2 =  
    end
    C = [1 -0.3014*exp(-j*2*pi*f) 0.3988*exp(-j*2*pi*f*2)];    
    D = (abs(sum(C)))^2;
    i = i+1; 
    
    Sd(i) = B/D;
end

Pdf = Deltaf*Sd;
Pd  = sum(Pdf);

% Estimação da DEP da seria simulada
se =dlmread('fGn_(4096)pts_alpha(0.8)_NrOfVM(1)_TopLev(9)_realiz3_filtroIIR2.txt');
[Pxx,w] = pwelch(se,[],[],(2*N));

Pxx = Pxx(2:N+1);

Pf = Deltaf*Pxx;
P  = sum(Pf);

% Calculo da DEP do modelo T-FARIMA  
ar_farima = [0.3014; -0.3988]; % sinais dos coeficientes de acordo com a convenção do S-PLUS = -1*Matlab

% Determina Polinomio "Vphi(B)=phi(B)D(B)" de acordo com a Eq. (16.26), pg. 475, Livro "Análise de
% Séries Temporais", Morettin e Toloi, onde D(B) = (1-B)^d = 1 -d(1)B -d(2)B^2 - ...

ar_farima = -1*ar_farima; % lembrar que os sinais dos coeficientes ajustados pelo S-PLUS são invertidos, por convenção
ar_farima = [1; ar_farima];

L=50;
M = length(ar_farima);
K = L-(M-1);

% calcula d(j)
for j = 0:K 
    D(j+1) = -gamma(j-d)/(gamma(j+1)*gamma(-d));
end

D = -1*D; % "D" é um vetor linha

Vphi_B = conv(ar_farima,D');
Vphi_B2 = -1*Vphi_B(2:length(Vphi_B)); 

Btfarima = 1;
Atfarima = Vphi_B;
F=(0.5/N):(0.5/N):0.5;
H = freqz(Btfarima,Atfarima,F,1);
St = abs(H).^2;

Ptf = Deltaf*St;
Pt  = sum(Ptf);

% Calculo da DEP do modelo AR estimado 
Bar = 1;
Aar = load('ar10h09.txt');
Aar = -1*Aar; % lembrar que os sinais dos coeficientes ajustados pelo S-PLUS são invertidos, por convenção
Aar = [1; Aar];
H2 = freqz(Bar,Aar,F,1);
Sar = abs(H2).^2;

Par = Deltaf*Sar;
Ptar  = sum(Par);

R=Pt/Pd;
R1=Pt/P;
R2=Pt/Ptar;

f=(0.5/N):(0.5/N):0.5;

figure %figura 1
plot(f,10*log10(St),'r')
xlabel('frequency')
ylabel('PSD (dB)')

hold on
pause

plot(f,10*log10(R*Sd),'g')

hold on
pause

plot(f,10*log10(R2*Sar),'g')

erro = 10*log10(St)-10*log10(R*Sd);

figure %figura 2
plot(f,10*log10(St)-10*log10(R*Sd))
xlabel('frequency')
ylabel('dB')

figure %figura 3
loglog(f,R*Sd,'g')
xlabel('frequency')
ylabel('PSD')

hold on
pause
loglog(f,St,'r')

hold on
pause
loglog(f,R1*Pxx,'m')

%hold on
%pause
%loglog(f,R2*Sar,'b')



