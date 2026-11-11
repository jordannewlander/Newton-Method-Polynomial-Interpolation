%Script for plotting the tenth degree (n = 10) interpolating polynomial for the Runge Function
%f(x) = 1/(1 + x^2) by both newtons and chebyshev method
clear All all

n_one = 10;
%we will use the table that we caluclated with the discrete whole numbers
%from -5 -> 5
x_one = -5 : 1 : 5;
y_one = [1/26 1/17 1/10 1/5 1/2 1 1/2 1/5 1/10 1/17 1/26];
a_one = Coef(n_one,x_one,y_one);

% point = -5 -> 5;
% Eval(n,x,a,point) to have more discrete values

xval = [-5:0.5:5];
for i = 1:length(xval)
    yval_one(i) = Eval(n_one,x_one,a_one,xval(i));
end

figure(1)
plot(x_one,y_one,'rx')
hold on
plot(xval,yval_one,'bo')
hold on
%am getting some weird behavoir on the ends of the domain, looks like they
%are diverging to pos infinity instead of converging to the x axis

%We will instead try doing it with Chebyshev's interpolation method, where
%we compute the xi values by using the formula 1/2(a+b) + 1/2(b-a)cos(i/n *
%pi) for i = 0 : 1 : n, a = -5, b = 5 and then compute coefficiants

for i = 0 : 10 
    x_two(i+1) = 1/2 * (-5 + 5) + 1/2 * (5 + 5) * cos(i/n_one * pi);
end
%now compute y values for this polynomial + coefficiants (coef)
y_two = 1./(1 + x_two.^2);
a_two = Coef(n_one, x_two, y_two);

for i = 1:length(xval)
    yval_two(i) = Eval(n_one,x_two,a_two,xval(i));
end

plot(x_two,y_two,'gx')
hold on
plot(xval,yval_two,'yo')
hold on
%this graph looks much better as it is converging to the x axis and not diverging anymore,
%why dont we also plot the original graph for comparison
xPlot = -5 : 0.1 : 5;
yPlot = 1 ./ (1 + xPlot.^2);
plot(xPlot, yPlot, 'w')

legend('data pnts newton','interpolation for newton', 'data pnts chebyshev', 'interpolation for chebyshev', 'original f(x)')
grid on
title('Interpolating Polynomial for Runge Function 1/(1 + x^2)')

%okay now plot diff g(x) values
%estimate is that as we have a higher degree polynomial (by increasing n), the approximmation for the
%polynomial will become increasingly more accurate
%n = 20
for i = 0 : 20 
    x_three(i+1) = 1/2 * (-5 + 5) + 1/2 * (5 + 5) * cos(i/20 * pi);
end
%now compute y values for this polynomial + coefficiants (coef)
y_three = 1./(1 + x_three.^2);
a_three = Coef(20, x_three, y_three);

for i = 1:length(xval)
    yval_three(i) = Eval(20,x_three,a_three,xval(i));
end

figure(2)
plot(x_three,y_three,'rx')
hold on
plot(xval,yval_three,'bo')
hold on

%n = 40
for i = 0 : 40 
    x_four(i+1) = 1/2 * (-5 + 5) + 1/2 * (5 + 5) * cos(i/40 * pi);
end
%now compute y values for this polynomial + coefficiants (coef)
y_four = 1./(1 + x_four.^2);
a_four = Coef(40, x_four, y_four);

for i = 1:length(xval)
    yval_four(i) = Eval(40,x_four,a_four,xval(i));
end

plot(x_four,y_four,'gx')
hold on
plot(xval,yval_four,'yo')
hold on

%n = 60
for i = 0 : 60 
    x_five(i+1) = 1/2 * (-5 + 5) + 1/2 * (5 + 5) * cos(i/60 * pi);
end
%now compute y values for this polynomial + coefficiants (coef)
y_five = 1./(1 + x_five.^2);
a_five = Coef(60, x_five, y_five);

for i = 1:length(xval)
    yval_five(i) = Eval(60,x_five,a_five,xval(i));
end

plot(x_five,y_five,'wx')
hold on
plot(xval,yval_five,'mo')

legend('data pnts g20','interpolation g20', 'data pnts g40', 'interpolation g40', 'data pnts g60', 'interpolation g60')
grid on
axis([-5 5 -0.5 2])
title('Interpolating Polynomial for Runge Function 1/(1 + x^2)')
%okay so I was wrong, as you add more data points the approximation does
%not get better, it stays the same until you have an approximation with too many nodes, 
%which also starts to diverge at the endpoints (g60 is starting to diverge
%at x = -5. you need to different methods of approximation 
%for different typres or rational polynomial functions (chebyshev vs
%newton, etc.), and you need to be mindful of the amount of nodes you are
%using to construct the polynomial
%the left endpoints for g60 is starting to diverge because the divided diff
%formula for xi = 5cos((i*pi)/n) (where n is the total number of nodes) is 
%ai+1(new) = (ai+1(old) - ai(old))/(xi+1 - xi+1-j) (found in coef function,
%used to find coefficients for each variable term in polynomial, derived in class)
%when two nodes are close together 