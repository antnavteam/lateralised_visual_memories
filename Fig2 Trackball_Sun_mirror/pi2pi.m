function th=pi2pib(th)
    th = mod(th + pi, 2*pi) - pi;
