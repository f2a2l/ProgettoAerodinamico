function [xd, yd, xv, yv] = NACA_5d (ddd, SS, n, varargin)
%
% NACA_5d.m
%
% Generazione dei punti sul dorso e sul ventre di un profilo alare NACA a 5
% cifre (five digits) a partire dalla linea media  linea_media_5d(x, ddd) e
% dallo spessore  spessore(x, SS)  del profilo simmetrico:
% 
% ddd : valore intero a tre cifre convenzionale    
%
% SS  : spessore massimo del profilo simmetrico   (due cifre)
%
% n   : numero di punti lungo la corda distribuiti uniformemente,
%       inclusi i punti del bordo di attacco e del bordo di uscita,  
%       in corrispondenza dei quali si calcolano i punti sul dorso
%       e sul ventre
%  
% Il parametro intero  ddd  con tre cifre puo' assumere solo
% i seguenti valori: 210, 220, 230, 240. 
%
% SS puo' essere fornito come NUMERO INTERO, 
% compreso fra 1 e 9,  in  %  della corda      
% oppure, se si preferisce, come VALORE REALE  < 1  
% frazione della corda
%             
%
% I. H. Abbot and A. E. von Doenhoff
% Theory of Wing Sections, pp 113
% Dover, New York, 1949, 1959
%


% Input
if SS >= 1
   s = SS/100; 
else  
   s = SS;
end

if (~isempty (varargin))
    spacing = varargin{1};
else
    spacing = 'constant';
end

xd = zeros(n+1);
yd = zeros(n+1);

xv = zeros(n+1);
yv = zeros(n+1);

m = n+1;

switch spacing
    
    case 'constant'
        
        for i = 1 : m
            
            x = (i-1)/(m-1);
            
            [y_lm, Dy_lm] = linea_media_5d(x, ddd);
            
            y_sp = spessore(x, s);
            
            theta = atan(Dy_lm);
            
            xd(i) = x  -  y_sp * sin(theta); 
            xv(i) = x  +  y_sp * sin(theta);
            
            yd(i) = y_lm  +  y_sp * cos(theta); 
            yv(i) = y_lm  -  y_sp * cos(theta);
            
        end
        
    case 'halfcos'
        
        n = n*2+1;
        m = n+1;
        i = [1:m];
        x = 1-0.5*(1+cos(((i-1)*pi)/n));
        x = x(1:ceil(length(x)/2))/0.5;
        
        for i = 1 : length(x)
            
            [y_lm, Dy_lm] = linea_media_5d(x(i), ddd);
            
            y_sp = spessore(x(i), s);
            
            theta = atan(Dy_lm);
            
            xd(i) = x(i)  -  y_sp * sin(theta); 
            xv(i) = x(i)  +  y_sp * sin(theta);
            
            yd(i) = y_lm  +  y_sp * cos(theta); 
            yv(i) = y_lm  -  y_sp * cos(theta);
            
        end
        
        if (xd(end,1) ~= 1)
            tmp = interp1 (xd(:,1), yd(:,1), 1, 'spline', 'extrap');
            xd(end+1,1) = 1;
            yd(end+1,1) = tmp;
        end
        
        if (xv(end,1) ~= 1)
            tmp = interp1 (xv(:,1), yv(:,1), 1, 'spline', 'extrap');
            xv(end+1,1) = 1;
            yv(end+1,1) = tmp;
        end
        
    case 'cos'
        
        for i = 1 : m
            
            x = 1-0.5*(1+cos(((i-1)*pi)/n));
            
            [y_lm, Dy_lm] = linea_media_5d(x, ddd);
            
            y_sp = spessore(x, s);
            
            theta = atan(Dy_lm);
            
            xd(i) = x  -  y_sp * sin(theta); 
            xv(i) = x  +  y_sp * sin(theta);
            
            yd(i) = y_lm  +  y_sp * cos(theta); 
            yv(i) = y_lm  -  y_sp * cos(theta);
            
        end
        
end

yd(end,1) = 0;
yv(end,1) = 0;

return

%==========================================================================
%==========================================================================
% linea_media_5d del profilo NACA a 5 cifre (five digits)

function [y, Dy] = linea_media_5d(x, ddd)

% Calcolo dell'ordinata della linea media  y(x)  e della
% sua derivata  dy(x)/dx  nel punto  x  lungo la corda
% del profilo NACA a 5 cifre con curvatura:  NACA-ddd__
%
% Le coordinate x e y sono adimensionali rispetto alla corda. 
% 
% Il parametro intero  ddd  con tre cifre puo' assumere solo
% i seguenti valori: 210, 220, 230, 240. 

% I. H. Abbot and A. E. von Doenhoff
% Theory of Wing Sections
% Dover, New York, 1949, 1959 

% Pagina 115-116

switch ddd

   case 210
     q = 0.0580;
     k = 361.4;
     
   case 220
     q = 0.1260;
     k = 51.64;
     
   case 230
     q = 0.2025;
     k = 15.957;
     
   case 240
     q = 0.2900;
     k = 6.643; 
     
   case 250
     q = 0.3910;
     k = 3.230;
     
   otherwise
     error ('Errore nella chiamata a linea_media_5d: la linea media richiesta non esiste.')
end


if x < q

    y = (k/6) * (x^3  -  3 * q * x^2  +  q^2 * (3-q) * x);

   Dy = (k/6) * (3 * x^2  -  6 * q * x  +  q^2 * (3-q)); 

else   

    y = (k/6) * q^3 * (1 - x);
  
   Dy = - (k/6) * q^3; 
   
end

return

%==========================================================================
%==========================================================================
% spessore di un profilo NACA (con 4 o 5 cifre)

function [y] = spessore(x, SS)

% Valore dello spessore  y  nel punto  x  lungo la corda
% per un determinato spessore massimo SS (due cifre)
%
% Le coordinate x e y sono adimensionali rispetto alla corda.
%
% Lo spessore massimo SS del profilo puo' essere fornito
% come NUMERO INTERO compreso fra 1 e 99  in  %  della corda
% oppure come VALORE REALE  < 1  come frazione della corda

% Ad esempio, entrambe le specificazioni  SS = 12  e  SS = 0.12 
% producono lo spessore del profilo simmetrico  NACA0012   

% I. H. Abbot and A. E. von Doenhoff
% Theory of Wing Sections
% Dover, New York, 1949, 1959 

% Pagina 113


if SS >= 1
   s = SS/100; 
else  
   s = SS;
end


y = 5 * s * (0.29690 * sqrt(x)  -  0.12600 * x  -  0.35160 * x^2 ...
                                +  0.28430 * x^3  -  0.10360 * x^4); %0.10150 * x^4);

return

%==========================================================================
%==========================================================================
% Derivata dello spessore di un profilo NACA (con 4 o 5 cifre)

function [Dy] = d_spessore_dx(x, SS)

% Valore della derivata dello spessore  y  nel punto  x  lungo 
% la corda per un determinato spessore massimo SS (due cifre)
%
% Le coordinate x e y sono adimensionali rispetto alla corda.
%
% Lo spessore massimo SS del profilo puo' essere fornito
% come NUMERO INTERO compreso fra 1 e 99  in  %  della corda
% oppure come VALORE REALE  < 1  come frazione della corda
% Ad esempio, entrambe le specificazioni  SS = 12  e  SS = 0.12 

% Ad esempio, entrambe le specificazioni  SS = 12  e  SS = 0.12 
% producono la derivata dello spessore del profilo simmetrico  
% NACA0012   

% I. H. Abbot and A. E. von Doenhoff
% Theory of Wing Sections
% Dover, New York, 1949, 1959 

% Pagina 113


if SS >= 1
   s = SS/100; 
else  
   s = SS;
end


Dy = 5 * s * (0.14845 / sqrt(x)  -  0.12600  -  0.7032 * x  ...
                                 +  0.8529 * x^2  - 0.4144*x^3);  %0.406 * x^3);


return