
%Weronika Hebda 
%wybrana gwiazda Wodnika - Sadalsuud
delta = (-5.571175); % deliknacja 
alfa = 21.52598 ; %Rektascensja

%wybrane przeze mnie miejsca 
%półkula N - warszawa
lambda_w=21;
fi_w=52;

%półkula S -buenos aires
lambda_b = -59;
fi_b = -34;

%okolice równika - Libreville
lambda_l = 9.5;
fi_l=0.25;

%przeliczanie na liczbe dni juliańskich
jd = julday(2001,1,24);


%przeliczenie daty na kalendarz juliański 
function jd = julday(y,m,d) 
    if m <= 2, y = y-1; m = m+12; end 
    jd = floor(365.25*(y+4716))+floor(30.6001*(m+1))+ d-1537.5;
end

%czas gwiazdowy Greenwich
function g = GMST(jd)
T = (jd-2451545)/36525;
g = 280.46061837 + 360.98564736629*(jd-2451545)+0.000387933*T.^2-T.^3/38710000;
g = mod(g,360); 
end

%Przeliczenie czasu słonecznego UT na czas gwiazdowy S 
%oraz obliczenie kąta godzinnego 
function [t] = katgodz(y,m,d,h,lambda,alfa) 
    jd = julday(y,m,d); %dni 
    g = GMST(jd); %stopnie 
    UT1 = h*1.002737909350795; %godziny 

    %obliczenie czasu gwiazdowego(w stopniach) 
    S = UT1*15 + lambda + g; 

    %obliczenie kąta godzinnego(w stopniach) 
    t = S - alfa*15; 
end   

