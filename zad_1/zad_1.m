%% Weronika Hebda

clear; 
a=6378137; 
e2=0.00669437999013; 

%wsp samolotu 
macierzDane= load('danelotu2.txt'); 
phi=macierzDane(:,1);  
lambda=macierzDane(:,2); 
h=macierzDane(:,3); 
rows = size(macierzDane, 1);

%wsp lotniska 
phiB=50.033333;
lambdaB=8.570556;
hB=111;

%wspolrzedne xyz samolotu
[xs,ys,zs] = geo2xyz(phi,lambda,h,a,e2);
%wspolrzedne xyz lotniska
[xl,yl,zl] = geo2xyz(phiB,lambdaB,hB,a,e2);

%xyz to neu
for i = 1:rows
    [n(i),e(i),u(i)] = xyz2neu(phiB,lambdaB,xs(i),ys(i),zs(i),xl,yl,zl);
end

%odleglosc skosna
s = sqrt(n.^2+e.^2+u.^2);

%odleglosc zenitalna
z = acosd(u./s);

%okreslanie momentu zajscia pod horyzont
for i = 1:rows
    idx = find(z >= 90);
    break
end

%azymuty
azymut = atand(e./n);
for i = 1:rows
    if n(i) > 0 && e(i) < 0
        azymut(i) = azymut(i) + 360;
    elseif n(i) < 0 && e(i) > 0
        azymut(i) = azymut(i) + 180;
    elseif n(i) < 0 && e(i) < 0
        azymut(i) = azymut(i) + 180;
    elseif azymut(i) > 360
        azymut(i) = azymut(i) - 360;
    end
end

%rysowanie trasy we współrzędnych geodezyjnych
geoscatter(phi,lambda,5,'.r');

%we wspolrzednych nue 
figure;
%wyswietlanie punktu za ktorym zachodzi za horyzontem
plot3(n,e,u, 'o-','MarkerFaceColor','red','MarkerEdgeColor','red','MarkerIndices',idx(1,1))
title('Trasa lotu w układzie horyzontalnym');
xlabel('n');
ylabel('e');
zlabel('u');
grid on

%we wspolrzednych xyz
figure;
plot3(xs,ys,zs); 
title('Trasa lotu w układzie kartezjańskim');
xlabel('x');
ylabel('y');
zlabel('z');
grid on

% geo to xyz
function [x,y,z] = geo2xyz(fi, lam, h, a, e2)
    N = a./sqrt(1-e2.*sind(fi).^2);
    x = (N+h).* cosd(fi) .* cosd(lam);
    y = (N + h) .* cosd(fi) .* sind(lam);
    z = (N.*(1-e2)+h).*sind(fi);
end

%xyz to neu
function [n,e,u] = xyz2neu(fi,lam,xA,yA,zA,xB,yB,zB)
    %macierz obrodu z wykorzystaniem wspolrzednych lotniska
    T = [-sind(fi).*cosd(lam) -sind(lam) cosd(fi).*cosd(lam);
     -sind(fi).*sind(lam) cosd(lam) cosd(fi).*sind(lam);
      cosd(fi) 0 sind(fi)];
    T=T';
    
    %wspl samolotu - wspl lotniska
    D = [xA - xB;yA - yB;zA - zB];

    neu = T*D;
    n = neu(1);
    e = neu(2);
    u = neu(3);
end
