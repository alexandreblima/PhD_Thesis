{\rtf1\ansi\ansicpg1252\deff0\deflang1033{\fonttbl{\f0\fnil\fcharset0 Courier New;}{\f1\fswiss\fcharset0 Arial;}}
{\*\generator Msftedit 5.41.15.1507;}\viewkind4\uc1\pard\f0\fs20 [Data, archive_name] = Model(TopLev, NrOfVM, h, g, model, varargin\{:\});\par
function [Data,nome_arquivo] = Model(TopLev, NrOfVM, h, g, modelo, varargin)\par
\tab   \par
    if modelo ~= 2, \par
        Lm = varargin\{1\};\par
        varargin(1) = []; \par
        switch modelo\par
            case 1, \par
                nome_arquivo = ['fGn_(',int2str(Lm),')pts_alpha(',num2str(varargin\{1\}),')_\par
                NrOfVM(',int2str(NrOfVM),')_TopLev(',int2str(TopLev),')'] ;\par
            case 3, \par
                nome_arquivo =  ['fGn_(',int2str(Lm),')pts_alpha(',num2str(varargin\{1\}),')_\par
                phi(',num2str(varargin\{2\}),')_NrOfVM(',int2str(NrOfVM),')_\par
                TopLev(',int2str(TopLev),')'] ;                 \par
            case 9, \par
                 nome_arquivo = [];\par
        end\par
        \par
    else %2 = MWM:\par
        if NrOfVM ~= 1,\par
            error('Modelo 2 (MWM) exige que se escolha NrOfVM = 1 para usar a Wavelet de Haar');\par
        end\par
        Lm = 2^TopLev; % one root scale coefficient\par
        \par
        nome_arquivo = ['MWM_(',int2str(Lm),')pts_alpha(',num2str(varargin\{1\}),')_p(',\par
        num2str(varargin\{2\}),')_TopLev(',int2str(TopLev),')'] ;\par
    end\par
    \par
    Lf = 2*NrOfVM; % filter length\par
\tab kn = [0 Lm-1]; \tab\par
\tab kp = kn; \par
\tab Data\{0+(1)\}.kp = kn; \par
\tab kc = kn;  \par
\tab kd = kn; \par
\tab Data\{0+(1)\}.kd = kn;\par
\tab\par
\tab for j=1:TopLev,\par
        kc = [kd(1)-(Lf-1) kd(2)]; \par
        kd = fix (kc/2); \par
         \par
        kn = floor([0 kn(2)-(Lf-1)]/2); \par
        kp = floor(kp/2) - [(NrOfVM - 1) 0]; \par
        \par
        Data\{j+(1)\}.kd = kd;\par
        Data\{j+(1)\}.kp = kp; \par
        Data\{j+(1)\}.kn = kn;\par
\tab end\par
\tab\par
\tab Data\{TopLev+(1)\}.app  = appProcess(diff(kp)+1, modelo, varargin\{:\}); \par
    \par
    for j=TopLev:-1:1,\par
        [Data] = detProcess(j,TopLev,Data,modelo,varargin\{:\});\par
        Data\{j-1+(1)\}.app = Recon(j,NrOfVM,Data,h,g); \par
    end\par
\par
    \par
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\par
%internal:\par
function approxs = appProcess(n, modelo, varargin)\par
\tab   switch modelo\par
        case 1, \par
            approxs = 4+randn(1,n);\par
        case 2, approxs = 4+randn(1,n); \par
        case 3, approxs = zeros(1,n); \par
        case 9, \par
            vetor_coefs = varargin\{1\};\par
            L=varargin\{2\};\par
            TopLev = varargin\{3\};\par
            wname = varargin\{4\};\par
            coefsA = appcoef(vetor_coefs,L,wname,TopLev);\par
            if length(coefsA) == n ,\par
                approxs = coefsA; \par
            else\par
                warning('Numero de coeficientes da dwt nao bate com o \par
                necessario para reconstituir o sinal');\par
                approxs = coefsA(1:n);\par
            end\par
        otherwise, disp('Modelo nao implementado'); approxs = []; \par
    end\par
\par
            \par
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\par
%internal:\par
function [Data] = detProcess(level, TopLev, Data, modelo, varargin)\par
    if modelo ~= 9,\par
        alpha = varargin\{1\}; \par
        ctrvariance = 2^( (level - TopLev) * alpha/2 );\par
    end\par
    \par
    n = diff( Data\{level+(1)\}.kp ) + 1;\par
    \par
    switch modelo\par
        case 1, \par
            Data\{level+(1)\}.det = ctrvariance*randn(1,n); \par
        case 2, %MWM:\par
            p = varargin\{2\};\par
            if (level == TopLev) & (p == 0),\par
                p = (2^alpha - 1)/(2 - 2^alpha);\par
                Data\{0+(1)\}.p = p;\par
                for i = 1:TopLev,\par
                    p = (2*p+1)/2^alpha - 1;\par
                    Data\{i+(1)\}.p = p;\par
                end\par
            elseif (level == TopLev) & (p > 0),\par
                Data\{TopLev+(1)\}.p = p;\par
                for i = TopLev-1:-1:0,\par
                    p = (2^alpha*(p+1)-1)/2;\par
                    Data\{i+(1)\}.p = p;\par
                end\par
                p = Data\{TopLev+(1)\}.p;\par
            elseif level ~= TopLev,\par
                p = Data\{level+(1)\}.p;\par
            end\par
            R01 = betarnd(p,p,1,n);\par
            R = 2.*R01-1;\par
            Data\{level+(1)\}.det = R.*Data\{level+(1)\}.app;\par
        case 3, phi = varargin\{2\}; \par
            inject = ctrvariance*randn(1,n);\par
            Data\{level+(1)\}.det(1) = phi*0 + inject(1);\par
            for i=2:n,\par
                Data\{level+(1)\}.det(i) = phi*Data\{level+(1)\}.det(i-1) + inject(i);\par
            end\par
            clear inject;\par
        case 9, \par
            vetor_coefs = varargin\{1\};\par
            L = varargin\{2\};\par
            coefsD = detcoef(vetor_coefs,L,level); \par
            if length(coefsD) == n ,\par
                Data\{level+(1)\}.det = coefsD; \par
            else\par
                warning('Numero de coeficientes da dwt nao bate com o \par
                necessario para reconstituir o sinal');\par
                Data\{level+(1)\}.det = coefsD(1:n);\par
            end     \par
            \par
        otherwise, error(['Modelo ', int2str(modelo), ' nao implementado!!!']); \par
    end\par
\par
\f1\par
}
 