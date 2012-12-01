%Recibe como inputs el largo L de la barra, el intervalo de distancia entre puntos h, intervalo de tiempo k, el factor alfa y el numero
%maximo de iteraciones a considerar como condicion de termino para el problema.

function [result] = eqHeatFD(L, h, k, alfa, MAX_ITER)
%Valores iniciales eqHeatFD(1,0.01,0.00000001,1,5000)	
	result = [];
	m = L/h;

%Se genera la partición de m intervalos con m+1 valores, correspondientes a los X_i a evaluar.
	X = linspace(0,L,m+1);
	u_0 = zeros(1,m+1);
	u_j = [];

%Se evaluan los X_i en la funcion de la condicion inicial para generar el vector U_0 inicial.
%De este vector, por condición de borde, el primer y último valor deben ser cero, por lo que no se evaluan
%en la función.
	for i = 2:m,
		u_0(i) = (exp(2*X(i)))*sin(X(i))*cos(X(i));
	end

%Se crea la matriz tridiagonal A.
	A = zeros(m+1,m+1);
	lambda = ((alfa^2)*k)/(h^2);
	dif = 1-(2*lambda);
	
	for i = 1:m+1,
		if i == 1
			A(1,i) = dif;
			A(1,i+1) = lambda;
		elseif i == m+1
			A(i,i) = dif;
			A(i,i-1) = lambda;
		else
			A(i,i) = dif;
			A(i,i-1) = lambda;
			A(i,i+1) = lambda;
	end
	
	u_ant = u_0;

	%Se ejecuta lo que aparece en el pseudo-codigo desde la linea 5
	%hasta la linea 12.
	for iter = 1:MAX_ITER,
		u_j = A*u_ant';
%Se aplica la condición de borde donde el primer y ultimo elemento de u_j deben ser cero.
		u_j(1) = 0;
		u_j(m+1) = 0;

		u_j2 = u_j';
		
		diferencia = u_j' - u_ant;
		
		%Se agregan vectores u_j generados a una matriz de resultados, para ser usada al momento de graficar.
		result(:,iter) = u_ant;
		
		u_ant = u_j2;		
		
		%Condición de término del algoritmo.
		if max(abs(diferencia)) <= (1.1*10^(-4))
			result(:,iter) = u_ant;
			break
		end
	end
end