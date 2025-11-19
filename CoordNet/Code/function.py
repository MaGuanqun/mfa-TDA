import torch



def range(function_name='expotential'):
    if function_name == 'expotential':
        return [-1.55, 1.55, -0.8, 2.3]
    elif function_name == 'schwefel':
        b=10.5 * torch.pi * 10.5 * torch.pi
        return [-b,b,-b,b]
    elif function_name == 'quartic_potential_2':
        return [-2.0, 2.0, -2.0, 2.0, 0.0, 4.0]
    else:
        raise NotImplementedError(f"Function {function_name} not implemented.")

def sample_size(function_name='expotential'):
    if function_name in ['expotential','schwefel']:
        return (401, 401)
    elif function_name in ['quartic_potential_2']:
        return (100, 100, 100)
    elif function_name == 'vortex_street':
        return (50, 80, 100)
    elif function_name == 'vortex_street_3d':
        return (640, 80, 150)
    elif function_name == 'hurricane_isabel':
        return (500, 500, 100)
    elif function_name == 'boussinesq_3d':
        return (150,450,200)
    elif function_name == 'fluid':
        return (512,512,100)
    elif function_name == 'cylinder':
        return (296,72,100)
    elif function_name == 'cylinder2':
        return (175,72,100)
    else:
        raise NotImplementedError(f"Function {function_name} not implemented.")

def function_2D(x,y, function_name='expotential'):
    if function_name == 'expotential':
        return compute_expotential(x, y)
    elif function_name == 'schwefel':
        return compute_schwefel(x, y)
    else:
        raise NotImplementedError(f"Function {function_name} not implemented.")

def function_3D(x, y, t, function_name='quartic_potential_2'):
    if function_name == 'quartic_potential_2':
        return compute_quartic_potential_2(x, y, t)
    else:
        raise NotImplementedError(f"Function {function_name} not implemented.")

def compute_expotential(x, y):
    term1 = torch.exp(-(8 * (x + 0.4) ** 2 + 4 * y ** 2))
    term2 = torch.exp(-8 * (x - 0.5) ** 2 - 4 * y ** 2)
    term3 = torch.exp(-8 * x ** 2 - 4 * (y - 0.77) ** 2)
    term4 = torch.exp(-8 * x ** 2 - 4 * (y - 1.5) ** 2)
    term5 = 0.2 * torch.exp(-0.3 * x ** 2 - 0.3 * (y - 0.5) ** 2)
    return term1 + term2 + term3 + term4 + term5

def compute_schwefel(x, y):
    """
    Schwefel function for 2D input (x, y).
    """
    a = 418.9829
    d = 2  # two dimensions: x, y

    term1 = x * torch.sin(torch.sqrt(torch.abs(x)))
    term2 = y * torch.sin(torch.sqrt(torch.abs(y)))

    result = 0.5 * (a * d - (term1 + term2))
    return result

def compute_quartic_potential_2(x, y, t):
    """
    Quartic potential function for 2D input (x, y) with time dependency.
    """
    return 0.25 * x**4 + 0.5 * (1 - t) * x**2 + 0.25 * y**4 + 0.5 * torch.cos(t) * y**2