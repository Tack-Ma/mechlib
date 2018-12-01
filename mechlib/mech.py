import numpy as np
import scipy.constants

def sec_data(**param):
    try:
        if param['sectype'] == 'rectangle':
            b = param['b']
            h = param['h']
            return b*h, b*h**3/12, b*h**2/6
        elif param['sectype'] == 'circle':
            d = param['d']
            return np.pi*d**2/4, np.pi*d**4/64, np.pi*d**3/32
        elif param['sectype'] == 'ring':
            d_o = param['d_o']
            d_i = param['d_i']
            return np.pi*(d_o**2-d_i**2)/4, np.pi*(d_o**4-d_i**4)/64, np.pi/32*(d_o**4-d_i**4)/d_o
    except:
        return np.nan,np.nan,np.nan

def beams_deflection(**param):
    try:
        if param['beamstype'] == '3.1':
            load = param['p']
            length = param['l']
            elasticity = param['E']
            I_y = param['I_y']
            x = param['x']
            omega = load*length**3/(6*elasticity*I_y)*(3-x/length)*x**2/length**2
            theta = load*length**2/(2*elasticity*I_y)*(2-x/length)*x/length
            return omega, theta
        elif param['beamstype'] == '3.2':
            load = param['p']
            length = param['l']
            elasticity = param['E']
            I_y = param['I_y']
            x = param['x']
            a = param['a']
            if x<=a:
                omega = load*a**3/(6*elasticity*I_y)*(3-x/a)*x**2/a**2
                theta = load*a**2/(2*elasticity*I_y)*(2-x/a)*x/a
            else:
                omega = load*a**3/(6*elasticity*I_y)*(2+3*(x/a-1))
                theta = load*length**2/(2*elasticity*I_y)*(2-x/length)*x/length
            return omega, theta
        elif param['beamstype'] == '3.3':
            return np.nan, np.nan
        elif param['beamstype'] == '3.4':
            qload = param['q']
            length = param['l']
            elasticity = param['E']
            I_y = param['I_y']
            x = param['x']
            omega = qload * length**2 / (24*elasticity*I_y) * x**2 * (6 - 4*x/length + x**2/length**2)
            theta = qload * length**2 / (6*elasticity*I_y) * x * (3 - 3*x/length + x**2/length**2)
            return omega, theta
    except:
        return np.nan, np.nan

def trq_to_kw( trq, rot,unit='Nm',):
    try:
        if unit == 'Nm':
            return 2 * np.pi * trq * rot / 60 / 1000
        elif unit == 'kgfm':
            return 2*np.pi*trq*scipy.constants.g * rot / 60 / 1000
        elif unit == 'kgfmm':
            return 2*np.pi*trq*scipy.constants.g /1000 * rot / 60 / 1000
    except:
        return np.nan

def kw_to_trq(p, rot, unit='Nm'):
    try:
        if unit == 'Nm':
            return 1000 * p / (2*np.pi)* 60/rot
        elif unit == 'kgfm':
            return 1000*p/(2*np.pi*scipy.constants.g)*60/rot
        elif unit == 'kgfmm':
            return 1000000*p/(2*np.pi*scipy.constants.g)*60/rot
    except:
        return np.nan
