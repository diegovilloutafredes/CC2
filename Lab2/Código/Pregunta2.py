#################################################################################################
import scipy.special as sps
import scipy.integrate as spi
import math as m
#################################################################################################
def get_roots(n):
#Obtiene las raices del Polinomio de Legendre de grado 'n'.  
    return sps.orthogonal.p_roots(n)[0]
#################################################################################################
def get_coefs(n):        
#Obtiene los coeficientes Ci asociados a cada una de las raices del Polinomio de Legendre de grado 'n'.
    return sps.orthogonal.p_roots(n)[1]
#################################################################################################
def erf_teo(y, n):
#Se comprueba que 'n' este dentro del rango pedido.
    if 4 <= n <= 7:

#Se obtienen las raices y los coeficientes necesarios para la cuadratura.
        x_i = get_roots(n)
        c_i = get_coefs(n)    

#La cuadratura de Gauss solo es valida en los limites de integracion [-1, 1], por lo que
#se deben cambiar los limites de la integral erf(y) a estos. Luego de realizar el cambio
#se obtiene un factor constante denominado 'factor_1'.        
        factor_1 = y/m.sqrt(m.pi)
        
        res = 0

#Es en este ciclo donde se realiza la cuadratura como tal, se hace una sumatoria de la
#multiplicacion de cada coeficiente Ci por la funcion evaluada en la raiz del Polinomio
#de Legendre correspondiente. En este caso, como se debieron mover los limites de integracion,
#la funcion no se evalua directamente en las raices y ahora se evalua en el valor contenido
#por la variable 'x'.
        for i in range(n):        
            x = (y*x_i[i] + y)/2        
            res += c_i[i]*m.exp(-x**2)

#Finalmente se multiplica el resultado obtenido de la sumatoria por el factor constante explicado
#anteriormente.
        return res*factor_1
    else: 
        print 'n value out of range'
#################################################################################################
def erf_real(y):
#Recibe como parametro el limite superior de la integral.

#'factor' representa la parte constante de la integral, y 'func' es la funcion a integrar.
    factor = 2/m.sqrt(m.pi)
    func = lambda x: m.exp(-x**2)

#Se utiliza scipy.integrate.quad() para obtener el valor de la integral definida entre [0, y].
    return factor*spi.quad(func,0,y)[0]
#################################################################################################
def err_teo_real(y,n):
#Retorna el error relativo entre el valor teorico y el real de erf(y) segun la formula explicada
#en el laboratorio.
    return abs(erf_teo(y, n) - erf_real(y))/erf_real(y)
#################################################################################################
#Funciones que solo imprimen por pantalla resultados de manera conveniente para
#ser usados en latex.
#################################################################################################
def print_roots():
    
    for i in range(4,8):
        roots = get_roots(i)
        
        if i == 4:
            print '$',i,'$ & $',roots[0],'$ & $',roots[1],'$ & $',roots[2],'$ & $',roots[3],'$ & $ $ & $ $ & $ $ \\\ '
        if i == 5:
            print '$',i,'$ & $',roots[0],'$ & $',roots[1],'$ & $',roots[2],'$ & $',roots[3],'$ & $',roots[4],'$ & $ $ & $ $ \\\ '
        if i == 6:
            print '$',i,'$ & $',roots[0],'$ & $',roots[1],'$ & $',roots[2],'$ & $',roots[3],'$ & $',roots[4],'$ & $',roots[5],'$ & $ $ \\\ '
        if i == 7:
            print '$',i,'$ & $',roots[0],'$ & $',roots[1],'$ & $',roots[2],'$ & $',roots[3],'$ & $',roots[4],'$ & $',roots[5],'$ & $',roots[6],'$ \\\ '

def print_coefs_erf():
    
    for i in range(4,8):        
        c_i = get_coefs(i)    
            
        if i == 4:
            print '$',i,'$ & $',c_i[0],'$ & $',c_i[1],'$ & $',c_i[2],'$ & $',c_i[3],'$ & $ $ & $ $ & $ $ & $',erf_teo(10, i),'$ \\\ '
        if i == 5:
            print '$',i,'$ & $',c_i[0],'$ & $',c_i[1],'$ & $',c_i[2],'$ & $',c_i[3],'$ & $',c_i[4],'$ & $ $ & $ $ & $',erf_teo(10, i),'$ \\\ '
        if i == 6:
            print '$',i,'$ & $',c_i[0],'$ & $',c_i[1],'$ & $',c_i[2],'$ & $',c_i[3],'$ & $',c_i[4],'$ & $',c_i[5],'$ & $ $ & $',erf_teo(10, i),'$ \\\ '
        if i == 7:
            print '$',i,'$ & $',c_i[0],'$ & $',c_i[1],'$ & $',c_i[2],'$ & $',c_i[3],'$ & $',c_i[4],'$ & $',c_i[5],'$ & $',c_i[6],'$ & $',erf_teo(10, i),'$ \\\ '

def print_error():
    
    for i in range(4, 8):
        print '$',i,'$ & $',erf_teo(10, i),'$ & $',erf_real(10),'$ & $',err_teo_real(10, i),'$ \\\ '

#################################################################################################
if __name__ == "__main__":
    
    for n in range(4,8):
        print erf_teo(10,n)