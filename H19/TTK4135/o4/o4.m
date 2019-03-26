clear all;
v = 0:0.01:6;
[X,Y] = meshgrid(v);
Z = 0.4*X.^2-3*X+0.2*Y.^2-2*Y;
conditions = (2*X+Y < 8) & (X+3*Y < 15);
cond = zeros(length(v)); % Initialize
cond(conditions) = NaN;
surface = surf(X, Y, cond);
surface.EdgeColor = 'g';
view(0,90)
hold on;
contour(X,Y,Z)
plot3(2.25,3.5,1,'x');
plot3(2.4,3.2,1,'o');
plot3(0,0,1,'o');

