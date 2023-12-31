{\rtf1\ansi\ansicpg1252\deff0\deflang1033{\fonttbl{\f0\fnil\fcharset0 Courier New;}{\f1\fswiss\fcharset0 Arial;}}
{\*\generator Msftedit 5.41.15.1507;}\viewkind4\uc1\pard\f0\fs20 function approxs = Recon(level, NrOfVM, Data, h, g) % h will be high-pass filter coefficients.\par
    Lf = length(h);\par
    \par
    kp = Data\{level + (1)\}.kp; \par
    \par
    % Eqs. (7.11), (7.12), (7.13) and (7.14) of Backar M.Sc thesis (see  p.49):\par
    %%%%%%%%%%%%%%%%%%\par
    \par
    ku = 2*kp; \par
    \par
    kc = ku + [0 Lf-1]; \par
    \par
    %%%%%%%%%%%%%%%%%%%\par
    \par
    kp = Data\{level-1 + (1)\}.kp; \par
    \par
    appro =  Data\{level + (1)\}.app;\par
    detail =  Data\{level + (1)\}.det;\par
    \par
    appup = dyadup(appro,0);\par
    detup = dyadup(detail,0);\par
    \par
    appconv = conv(appup,g); \par
    detconv = conv(detup,h); \par
    \par
    indices = [kp(1)-kc(1)+(1):kp(2)-kc(1)+(1)];\par
    \par
    appro = appconv(indices) + detconv(indices);\par
    \par
    approxs = appro;\par
\f1\par
}
 