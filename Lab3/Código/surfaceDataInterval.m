function [] = surfaceDataInterval(L, MAX_ITER)
%L = 0.5, MAX_ITER = 20
	result = eqHeatFD(L, 0.01, 0.00000001, 1, MAX_ITER);
	plot(result')
end