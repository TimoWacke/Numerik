'''
    Newton Interpolation
    and inheriting from that Hermite Interpolation and Tschebyscheff Interpolation

    Cubic Spline Interpolation is also implemented
'''


from scipy.interpolate import CubicSpline
import numpy as np


# import classes from other files
from divided_differences import DividedDifferences
from divided_differences import SupportingPoint
from divided_differences import SupportPointList


class Newton:
    def __init__(self, spl: SupportPointList):
        self.spl = spl
        self.dd_calc = DividedDifferences(spl)
        self.dds = []

    def basisPolynom(self, n, x):
        # yields n times up to the basis with grade n
        assert n < len(self.spl.list)
        product = 1
        yield product
        for i in range(1, n+1):
            product *= (x - self.spl.get(i-1).x)
            yield product

    def interpolate(self, n=False, visualize=False):
        if not n:
            n = len(self.spl.list) - 1
        assert n < len(self.spl.list) and n >= 0
        local_dds = [] # first grade
        for i in range(0, n):  # 0, 1, ..., n-1, next line n-th grade
            local_dds.append(self.dd_calc.dds_nth_grade(i))  # not visualizing the dds   
        local_dds.append(self.dd_calc.dds_nth_grade(n, visualize=visualize))
        self.dds = local_dds
        # using generator function for the basis polynom
        # and using saved list of dds

        return lambda x: sum([dd * product for dd, product in zip(self.dds, self.basisPolynom(n, x))])

'''
    Class that can interpolate a given function f using Hermite Interpolation
    and the DividedDifference class
    Inherits from Newton as most of the process is the same
'''

class Hermite(Newton):
    def __init__(self,f):
        self.f = f
        # the rest can be achieved with Newton's init
        super().__init__(SupportPointList())
    
    def generateDataPoints(self, n, m, a=-1, b=1):
        x_list = np.linspace(a, b, n)
        for x in x_list:
            self.addDataPoints(x, m)

    def addDataPoints(self, x, m=1): # use the first m-1 derivatives as data point, m=1 means using just f
        f_at_x = self.f(x)
        derivative_suports = list(FunctionOverX(f).n_th_derivative_at(x, m))
        k_base = len(self.spl.list)
        self.dd_calc.expand_size_to(k_base + m)
        for j in range(m):
            self.spl.add(SupportingPoint(x, f_at_x))
            for k in range(j, m):
                self.dd_calc._k_js[k_base + k][j] = derivative_suports[j]                 


'''
    Class that can interpolate a given function f using Tschebyscheff Interpolation
    and the DividedDifference class
    Inherits from Newton as most of the process is the same

'''        

class Tschebyscheff(Newton):
    def __init__(self, f):
        self.f = f
        # the rest can be achieved with Newton's init
        super().__init__(SupportPointList())

    def generateDataPoints(self, n):
        # generate n data points
        # x_i = cos((2i+1)/(2n+2) * pi)
        # i = 0, ..., n
        for i in range(n):
            x_i = np.cos((2*i+1)/(2*(n-1)+2) * np.pi)
            self.spl.add(SupportingPoint(x_i, self.f(x_i)))
            self.dd_calc.expand_size_to(i+1)
        
        


'''
    Wrapper class for scipy's CubicSpline based on a given function f
'''
class Spline:

    def __init__(self, f):
        self.f = f # function that we want to interpolate
        self.x_list = None
        self.y_list = None
        self.spline = None

    def generateDataPoints(self, n, a=-1, b=1):
        # generate n data points
        # x_i = a + i * (b-a) / n
        # i = 0, ..., n
        self.x_list = np.linspace(a, b, n)
        self.y_list = self.f(self.x_list) # runs over all x_i and calculates f(x_i)
        

    def interPolateWithCubicSpline(self):
        # natural cubic spline:
        # -> meaning
        # f''(x_0) = 0 and f''(x_n) = 0 
    
        self.spline = CubicSpline(self.x_list, self.y_list, bc_type='natural')

