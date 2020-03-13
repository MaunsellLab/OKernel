box = zeros(1,100);
box(55:65) = 1;

delta = zeros(1,10);
delta(5) = 1;

plot(conv(box,delta,'same'));

