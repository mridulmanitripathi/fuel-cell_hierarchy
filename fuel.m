function dydx = fuel(x,y)
w = 100;
Da = 61.75;
dydx = [y(2); i*w*y(1)/Da];
end