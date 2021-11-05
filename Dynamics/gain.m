
function g = gain(x,theta)

g = (x<theta).*0 + (x>=theta).*tanh((x-theta)*2/3);
