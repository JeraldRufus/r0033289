filename = 'trial_40degree (1)_1603.Direct simulation.1.txt';
x=readtable(filename);
p=removevars(x,{'Var8','Var10'});
y=table2array(p);
new1=rand(size(y,1),1);
new2=rand(size(y,1),1);

y=[y new1];
y=[y new2];

for i=1:1:size(y,1)
    y(i,9)= (atan(y(i,7)/y(i,5)))*180/3.14;
    y(i,10)= (atan(y(i,6)/y(i,5)))*180/3.14;
end

%[g,h]=meshgrid(y(:,9),y(:,10));
%z=g+h;
f=y(:,8);g=y(:,9);h=y(:,10);
%f=ones(1000,1000);
%surf(g,h,f)

%plot3(g,h,f,'o');
x = linspace(min(g)-10, max(g)+10,100);
y = linspace(min(h)-10,max(h)+10,numel(x)) ;
[XX, YY] = meshgrid(x,y);
%F = griddedInterpolant(g,h,f,'linear','none');
F = scatteredInterpolant(g,h,f,'linear','none');
ZZ = F(XX,YY);
mesh(XX,YY,ZZ);
%ZZ=griddata(g,h,f,XX,YY);
figure(2);
contour(XX,YY,ZZ,100);
%figure(3);
%[c,h]=contourm(XX,YY,ZZ,100);
%clegendm(c,h);
op=-40;
z=interp2(XX,YY,ZZ,op,YY(:,1));
figure(4);
plot(YY(:,1),z)




