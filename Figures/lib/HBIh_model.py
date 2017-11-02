# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:46:06 2016

@author: keshengXu
"""
import numpy as np
from scipy import integrate


"""
Aquí definimos una funcion que toma las variables a tiempo t
y devuelve las derivadas
"""
def HBIh(time, gsd = None, gh = None):
    """
    paras
    
    gsd: optional 
        Conductance of slow depolarizing
    
    gh : optional 
        
    
    Return
    
    Var_t : action potential 

    """
    
    def HyB(Var,t,tempF):
        [rrho,pphi] = tempF 
        [v,ar,asd,ca,ah]=Var
        ad = 1/(1+np.exp(-zd*(v-V0d)))
        isd = rrho*gsd*asd*(v - Ed)
        #Imemb=isd + rho*gd*ad*(v - Ed) + rho*(gr*ar + gsr*asr)*(v-Er) + gl*(v - El)
        Imemb = isd + rrho*gd*ad*(v - Ed) + rrho*(gr*ar + gsr*(ca**2)/(ca**2+0.4**2))*(v-Er) + rrho*gl*(v - El) \
                    + rrho*gh*ah*(v - Eh) 
        arinf = 1/(1+np.exp(-zr*(v-V0r)))
        asdinf = 1/(1+np.exp(-zsd*(v-V0sd)))
        ahinf = 1/(1+np.exp(-zh*(v-V0h)));
        
        return np.array([-Imemb,
                    pphi*(arinf - ar)/tr,
                    pphi*(asdinf - asd)/tsd,
                    pphi*(-eta*isd - kappa*ca)/tsr,
                    pphi*(ahinf-ah)/th])
    
    #Parámetros del modelo
    gd = 2.5; gr = 2.8; gsr = 0.28;
    gl = 0.06;
    V0d = -25; V0r = -25; zd = 0.25; zr = 0.25;tr = 2;
    V0sd = -40; zsd = 0.11; tsd = 10;
    eta = 0.014; kappa = 0.18; tsr = 35;
    V0h= -85; zh = -0.14; th=125;
    
    
    Ed = 50; Er = -90; El = -80; Eh = -30;
    
    if gsd is None and gh is None:
        gsd = 0.21
        gh = 0.4

    temp=36
    rho=1.3**((temp-25.)/10)
    phi = 3**((temp-25.)/10)
    
    #voltaje inicial e inicialización de variables
    v = -60
    #Luego calulamos el valor de las variables a ese voltaje
    ad = 1/(1+np.exp(-zd*(v-V0d)));
    ar = 1/(1+np.exp(-zr*(v-V0r)));
    asd = 1/(1+np.exp(-zsd*(v-V0sd)));
    ca = -eta*rho*gsd*asd*(v - Ed)/kappa;
    ah = 1/(1+np.exp(-zh*(v-V0h)));
    
    #Ahora viene la simulacion misma
    #Creamos un vector con los valores iniciales
    X=np.array([v,ar,asd,ca,ah])
    
    Var_t = integrate.odeint(HyB, X, time, args = ((rho,phi),))
    
    return Var_t







