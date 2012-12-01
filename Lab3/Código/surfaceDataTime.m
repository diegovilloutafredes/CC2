function [] = surfaceDataTime(MAX_ITER)
%Recibe el parametro MAX_ITER el cual segun el laboratorio debe ser 5000
%y corresponde al numero maximo de iteraciones del algoritmo antes de detenerse.
	mesh(eqHeatFD(1,0.01,0.00000001,1,MAX_ITER))
end