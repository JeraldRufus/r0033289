filename = 'trial_40degree (1)_1603.Direct simulation.1.txt';
x=readtable(filename);
p=removevars(x,{'Var8','Var10'});%,'Var11','Var12','Var13','Var14'});
y=table2array(p);
new1=rand(size(y,1),1);
new2=rand(size(y,1),1);

y=[y new1];
y=[y new2];


for i=1:1:size(y,1)
    y(i,9)= (atan(y(i,7)/y(i,5)))*180/3.14;
    y(i,10)= (atan(y(i,6)/y(i,5)))*180/3.14;
end

scatter(y(:,9),y(:,10))
xlim([-60 60]);
ylim([-60 60]);
[mean(y(:,9));mean(y(:,10))]