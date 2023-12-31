{\rtf1\ansi\ansicpg1252\deff0\deflang1033{\fonttbl{\f0\fnil\fcharset0 Courier New;}{\f1\fswiss\fcharset0 Arial;}}
{\*\generator Msftedit 5.41.15.1507;}\viewkind4\uc1\pard\f0\fs20 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\par
%                                                       %\par
%   Generate.m                                          %\par
%                                                       %\par
%   Fernando L. de Mello                                %\par
%   Alexandre B. de Lima                                %\par
%                                                       %\par
%   10/2006                                             %\par
%   12/2007                                             %\par
%                                                       %\par
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\par
%   References:\par
%    1) J.-A. B\'e4ckar, "A Framework for Implementing Fractal Traffic Models \par
%    in Real Time", M.Sc. Thesis, SERC, Melbourne,2000.\par
%    2) F. L. de Mello, "Estudo e Implementa\'e7\'e3o de um Gerador de Tr\'e1fego \par
%    com Depend\'eancia de Longa Dura\'e7\'e3o", M.Sc. Thesis, Univ. of Sao Paulo, 2006.   \par
%    3) F. L. de Mello, A. B. de Lima, M. Lipas, J. R. A. Amazonas, \par
%    "Gera\'e7\'e3o de S\'e9ries Auto-Similares Gaussianas via Wavelets para Uso em \par
%    Simula\'e7\'f5es de Tr\'e1fego", IEEE Latin America Transactions, Mar\'e7o, 2007.\par
%    4) A. B. de Lima, F. L. de Mello, M. Lipas, J. R. A. Amazonas,\par
%    "A Generator of Teletraffic with Long and Short-Range Dependence",\par
%    12th CAMAD Workshop (part of the 18th IEEE PIMRC07), Greece, 2007.\par
%    5) D. B. Percival and A. T. Walden, "Wavelet Methods for Time Series \par
%    Analysis", Cambridge University Press, 2000.\par
%    6) R. H. Riedi, M. S. Crouse, V. J. Ribeiro, R. G. Baraniuk, \par
%    "A multifractal wavelet model with application to network traffic"\par
%    IEEE Transactions on Information Theory, 45(3), pp.992-1018, April, 1999.\par
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \par
%\par
%    'Generate.m' is a wavelet-based generator of self-similar sample paths of \par
%    Fractional Gaussian Noise (FGN) and Multifractal Wavelet Model (MWM). \par
%    It works in conjunction with functions 'Model.m' and 'Recon.m'\par
%\par
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\par
% Models:\par
% 1 - wavelet coeficients as white noise random processes (FGN synthesis)\par
% 2 - MWM\par
% 3 - wavelet coeficients as AR(1) random processes (FGN synthesis)\par
% 9 - Reconstrucao Analise Matlab\par
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\par
% Parameters:\par
% TopLev    % number of iterations (reconstruction using Pyramid algorithm)\par
% NrOfVM    % number of vanishing moments of the Daubechies wavelet function\par
% modelo    % stochastic model used \par
% Lm        % number of generated samples\par
% alpha     % scaling exponent; remember that H = (1+alpha)/2 and d = alpha -0.5 \par
%           % where 'H' is the Hurst parameter and 'd' is the fractional parameter\par
%           % of an ARFIMA(p,d,q) model\par
% p         % shape parameter of the beta probability density function \par
%           % In the MWM model, the variance of the multiplier R is given by\par
%           % var[R] = 1/(2p+1). \par
% phi       % AR(1) coefficient\par
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\par
%\par
% use as: [Data, nome_arquivo]=Generate(TopLev, NrOfVM, modelo, varargin)\par
%      or [Data]=Generate(TopLev, NrOfVM, 9, varargin) - modelo 9\par
% where: \par
% for 'modelo'= 1 (FGN, white noise wav. coef.), 'varargin' will be 2 parameters:  "Lm, alpha"\par
% for 'modelo'= 2 (MWM), 'varargin' sera dois parametros:    "alpha, p"  \par
% for 'modelo'= 3 (FGN, AR(1) wav. coef.), 'varargin' will be 3 parameters: "Lm, alpha, phi"\par
% for 'modelo'= 9, 'varargin' will be five parameters:  "Lm, vetor_coefs (C), L, TopLev, wname" \par
%\par
% Example of interactive use: generation of an FGN signal with 4096 samples\par
%                             alpha = 0.6 (H=[1+alpha]/2=0.8), \par
%                             from 1 root scale (approximation) coeficient\par
%\par
% >> [Data, archive_name] = Generate(12, 1, 1, 4096,0.6)\par
% \par
% wname =\par
% db1\par
% g =\par
%    0.7071    0.7071\par
% Lf =\par
%    2\par
% h =\par
%    0.7071   -0.7071\par
% Data = \par
% Columns 1 through 6\par
%    [1x1 struct]    [1x1 struct]    [1x1 struct]    [1x1 struct]    \par
%[1x1 struct]    [1x1 struct]\par
%  Columns 7 through 12\par
%    [1x1 struct]    [1x1 struct]    [1x1 struct]    [1x1 struct]\par
%[1x1 struct]    [1x1 struct]\par
%  Column 13\par
%    [1x1 struct]\par
% archive_name =\par
% fGn_(4096)pts_alpha(0.6)_NrOfVM(1)_TopLev(12)\par
%\par
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\par
\par
function [Data, archive_name] = Generate(TopLev, NrOfVM, model, varargin)\par
\par
N = int2str(NrOfVM); % N is a string\par
wname = 'db'; % must concatenate N obtain dbN\par
% Compute the corresponding scaling filter. \par
wname = strcat(wname,N)\par
% Compute the corresponding scaling filter. \par
g = dbwavf(wname); \par
g = g * sqrt(2)\par
Lf = length(g) % must be 2N\par
h = fliplr(g); h = (-1).^(0:Lf-1).*h \par
% Eq.(75b) of book "Wavelet Methods for Time Series Analysis"\f1\par
}
 