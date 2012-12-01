#################################################################################################
import numpy as np
import scipy as sp
import timeit as t
#################################################################################################
def generate_Data(n):
#Se crea el vector X con n datos equiespaciados entre si, que pertenecen al rango [-1,1].
    x = sp.linspace(-1,1,n)
    y = []    

#Se genera el vector Y con n datos, correspondientes a evaluar los valores de X en la funcion
#original.
    for i in range(n): 
        
        num = float(10*sp.log10(x[i]**2 + x[i] + 1))
        den = float(10*(x[i]**3) - 20*(x[i]**2) + x[i] - 2)
        y.append(float(num/den))
    
    return x, y
#################################################################################################
def get_Yreal(x):
#Se evalua la funcion original en el valor x recibido.
    num = float(10*sp.log10(x**2 + x + 1))
    den = float(10*(x**3) - 20*(x**2) + x - 2)
    return float(num/den)
#################################################################################################
def diff(x_int, x, y, n):
#Se crea un vector que contiene los mismos valores de Y, y se mantiene su primer valor ya que
#este corresponde al primer coeficiente del polinomio.
    coef = y
    n = n-1

#El metodo de las diferencias divididas se aplica a continuacion, donde se van restando y
#dividiendo los valores del vector 'coef' y de X hasta que todos los valores contenidos en
#'coef' corresponden a los coeficientes finales del polinomio interpolador.
    for i in range(1, n+1):
        for j in range(n, i-1, -1):
            coef[j] = (coef[j] - coef[j-1])/(x[j] - x[j-i])

#Finalmente se evalua el valor de x_int en el polinomio de la forma
#coef[0] + coef[1](x_int - x[0]) + coef[2](x_int - x[0])(x_int - x[1])...
    resultado = coef[0]
    factor = 1
    
    for i in range(0, n):
            
        factor *= (x_int - x[i])
        resultado += factor*coef[i+1]
    
    return resultado
#################################################################################################
def splines(x_int, x, y, n):
#Se crea una matriz de ceros, de dimensiones 'nxn' que contiene los terminos [1, 4, 1]
#a lo largo de su diagonal a excepcion de las casillas [0,0] y [n-1, n-1].
    A = sp.zeros((n,n))
    A[0,0] = 1
    A[n-1,n-1] = 1
    
    for i in range(1,n-1):
        A[i,i-1] = 1
        A[i,i] = 4
        A[i,i+1] = 1
    
#Se crea un vector de ceros, de la misma dimension que el vector Y, el cual tiene ceros
#en el primer y ultimo elemento, y el resto del vector posee elementos de la forma
#y[i-1] - 2*y[i] + y[i+1].
    Y = sp.zeros(n)
    
    for i in range(1,n-1):
        Y[i] = y[i-1] - 2*y[i] + y[i+1]

#Una vez que se tienen ambos elementos, A e Y, se procede a resolver el sistema de ecuaciones
#con la funcion linalg.solve(A,Y) de numpy, la que retorna finalmente un vector con todos los
#coeficientes C_j de los splines.
    C = np.linalg.solve(A,Y)

#Finalmente se debe evaluar el valor de x_int en el spline correspondiente, para encontrar
#el coeficiente C_j que le corresponde, se busca primero el indice del vector X del valor mas
#cercano a x_int y luego que se tiene este indice, se procede a encontrar los demas coeficientes
#A_j, B_j y D_j utilizando las formulas que aparecen en la pagina de la asignatura.
    index = min(range(n), key=lambda i: abs(x[i]-x_int))

#El siguiente IF comprueba que el indice del spline no se salga de los limites del arreglo
#pudiendo asi evaluar el punto x_int.
    if index == (n-1):
        index = index - 1
    
    h = float((x[n-1] - x[0])/n)
    x_j = x[index]  
    a_j = y[index]
    c_j = C[index]    
    b_j = (y[index+1]-a_j)/h - (h*(C[index+1]+2*c_j))/3
    d_j = (C[index+1] - c_j)/(3*h)
    
	resultado = a_j + b_j*(x_int - x_j) + c_j*(x_int - x_j)**2 + d_j*(x_int - x_j)**3
    
    return resultado
#################################################################################################
def inter_pol(x_int, n, pol):
#Se generan los vectores X e Y que contienen los datos para realizar la interpolacion.
    x,y = generate_Data(n)

#Luego se llama a las funciones interpoladores disponibles segun se haya escogido en un principio.
    if pol == 'diff':        
        return diff(x_int, x, y, n)
            
    if pol == 'spl':        
        return splines(x_int, x, y, n)
#################################################################################################
def benchmark():
    
    X = [-0.5,-0.25,0.0,0.25,0.5]
    interp_methods = ['diff','spl']
    
    for m in range(1,6):
        for j in interp_methods:
            for i in X:
                
                n = 2**m
                y_real = get_Yreal(i)
                y_int = inter_pol(i,n,j)
                
                if y_real != 0:
                    error_rel = abs((y_real - y_int)/y_real)
                else:
                    error_rel = abs(y_real - y_int)
                
                comp_time = t.timeit("inter_pol(x_int,n,pol)",setup='x_int='+str(i)+'; n='+str(n)+'; pol="'+j+'"; from __main__ import inter_pol',number=10)
                
                if y_real == 0:
                    print 'N:',n,' X:',i,' Y_k:',y_real,'            Y_int:',y_int,' pol:',j,' Error Relativo:',error_rel,' Tiempo de computo:',comp_time
                else:
                    print 'N:',n,' X:',i,' Y_k:',y_real,' Y_int:',y_int,' pol:',j,' Error Relativo:',error_rel,' Tiempo de computo:',comp_time
            print ''
#################################################################################################
if __name__ == "__main__":

    benchmark()