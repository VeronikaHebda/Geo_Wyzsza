
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

godz = [1:25];
%dla Warszawy%
[xw,yw,zw,wysw,azymutw] = oblicz(fi_w,lambda_w,alfa,delta);

%wykres wysokości od czasu
figure(1)
plot(godz,wysw)
grid on
title('Wykres zależności wysokości od czasu dla Warszawy ')
xlabel('czas [h]');
ylabel('Wysokość gwiazdy nad horyzontem [°]');
xticks([0:24])

%wykres azymutu od czasu 
figure(2)
plot(godz,azymutw)
grid on
title('Wykres zależności azymutu od czasu dla Warszawy')
xlabel('czas [h]');
ylabel('azymut [°]');
xticks([0:24])


%rysowanie sfery
figure(3)
[x,y,z] = sphere(50); 
z(z<0) = 0;
S = surf(x,y,z); 
S.EdgeColor = 'blue';
S.FaceAlpha = 0.4;
light;
lighting gouraud;
axis equal;
hold on 
scatter3(xw,yw,zw, 'r','filled');
title('Pozorna wędrówka gwiazdy Sadalsuud dla Warszawy 24.01.2001')
xlabel('x');
ylabel('y');
zlabel('z');



%funkcja na obliczanie wspolrzednych horyzontalnych
function [x,y,z,wys,azymut] = oblicz(fi,lam,alfa,delta) 
    for i = 1:25
        hour(i) = katgodz(2001,01,24,i-1,lam,alfa);
        zenit(i) = acosd(sind(fi)*sind(delta)+cosd(fi)*cosd(delta)*cosd(hour(i)));
        a =(-cosd(delta)*sind(hour(i)));
        b =(cosd(fi)*sind(delta)-sind(fi)*cosd(delta)*cosd(hour(i)));
        azymut(i) = atand(a/b)
    
        wys(i) = 90 - zenit(i);
    end
    for i = 1:24
        x(i) = sind(zenit(i))*cosd(azymut(i)); 
        y(i) = sind(zenit(i))*sind(azymut(i)); 
        z(i) = cosd(zenit(i)); 
    end
end


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

