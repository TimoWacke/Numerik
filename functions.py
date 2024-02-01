'''
    class that can run backpropagation on a MathExpression to help us
    compute the derivative of a function at a given point
'''

import numpy as np
import matplotlib.pyplot as plt


class MathExpression:
    # value is the result of the math expression
    # children can be other MathExpressions or numbers
    # operator can be used for printing / representation
    def __init__(self, value, children=(), operator=''):
        self.value = value
        self.children = children
        self.operator = operator
        '''
        grad will be set by backpropagation as the partial derivative,
        of this expression relative to one of its later parents
        '''
        self.grad = 0 
        self._backward = lambda: None 

    def __add__(self, other):
        if not isinstance(other, MathExpression):
            other = MathExpression(other, (), 'const')
        out = MathExpression(self.value + other.value, (self, other), '+')
        
        def backward():
            self.grad += out.grad
            other.grad += out.grad
        out._backward = backward

        return out

    def __mul__(self, other):
        if not isinstance(other, MathExpression):
            other = MathExpression(other, (), 'const')
        out = MathExpression(self.value * other.value, (self, other), '*')

        def backward():
            self.grad += other * out.grad
            other.grad += self * out.grad

        out._backward = backward
        return out    

    def __pow__(self, other):
        # other needs to be constant
        assert isinstance(other, (int, float))
        other = MathExpression(other, (), 'const')
        out = MathExpression(self.value ** other.value, (self, other), '**')

        def backward():
            self.grad += other * (self ** (other.value - 1)) * out.grad
        out._backward = backward
        return out
    
    def sin(self):
        out = MathExpression(np.sin(self.value), (self,), 'sin')

        def backward():
            self.grad += out.grad * self.cos()
        out._backward = backward
        return out
    
    def cos(self):
        out = MathExpression(np.cos(self.value), (self,), 'cos')

        def backward():
            self.grad += out.grad * -self.sin()
        out._backward = backward
        return out
    
    def __neg__(self):
        return self * -1
    
    def __sub__(self, other):
        return self + (-other)
    
    def __radd__(self, other):
        return self + other
    
    def __rsub__(self, other):
        return -self + other
    
    def __rmul__(self, other):
        return self * other

    def exp(self):
        out = MathExpression(np.exp(self.value), (self,), 'exp')

        def backward():
            self.grad += out.grad * self.exp()
        out._backward = backward
        return out
    
    def backward(self):
         # topological order all of the children in the graph
        topo = []
        visited = set()
        def build_topo(v):
            if v not in visited:
                visited.add(v)
                for child in v.children:
                    build_topo(child)
                topo.append(v)
        build_topo(self)

        # go one variable at a time and apply the chain rule to get its gradient
        self.grad = MathExpression(1, (), 'const')
        for v in reversed(topo):
            v._backward()
  
    
    def reset_grad(self):
        self.grad = 0
        for child in self.children:
            child.reset_grad()
    
    def __repr__(self):
        if self.operator == 'const':
            return f'{self.value}'
        elif self.operator.startswith('var='):
            return f'Math({self.operator.split("=")[1]})={self.value}'
        elif len(self.children) == 1:
            return f'Math({self.operator}{self.children[0]})={self.value}'
        elif len(self.children) == 2:
            return f'Math({self.children[0]} {self.operator} {self.children[1]})={self.value}'
        else:
            return f'Math({self.value})'

'''
    class that enables derivation of a Functional at a given point
    using the MathExpression class which implements backpropagation
'''
class FunctionOverX:
    def __init__(self, function):
        self.function = function

    def __call__(self, x):
        x = MathExpression(x, (), 'var=x')
        return self.function(x).value
    
    '''
        calulates (d/dx)^n at self.f(x), arguments: x, n:int 
    '''
    def first_n_th_derivative_at(self, x, n=1):
        x = MathExpression(x, (), 'var=x')
        y = self.function(x)
        i_th_derivative = y
        yield i_th_derivative.value
        for _ in range(n):
            i_th_derivative.reset_grad()
            i_th_derivative.backward()
            i_th_derivative = x.grad
            yield i_th_derivative.value
    
    def n_th_derivative(self, n=1):
        return lambda x: list(self.first_n_th_derivative_at(x, n))[-1]

'''
    Class that makes it easy to plot multiple callables
    over a given interval with a given resolution

    Usage Example:
    plotter = PlotFunctions([('f', f), ('spline', spline.spline, x_supports, y_support)], spline.x_list, spline.y_list)

    x_supports and y_supports are optional and can be used to plot the support points of the spline
'''  

class PlotFunctions:
    def __init__(self, labled_functions):
        self.labled_functions = []
        for f in labled_functions:
            assert len(f) >= 2
            assert isinstance(f[0], str)
            assert callable(f[1])
            if len(f) == 2:
                self.labled_functions.append((f[0], f[1], None, None))
            else:
                assert len(f) == 4
                self.labled_functions.append((f[0], f[1], f[2], f[3]))


    
    def plot(self, a=-1, b=1, resolution=10):
        # resolution: drawn points per unit
        # xs length
        n = int(b - a) * resolution
        xs = np.linspace(a, b, n)
        # plot the function
        for label, function, xp, yp in self.labled_functions:
            # generate y values
            ys = []
            for x in xs:
                y = function(x)
                ys.append(y)
            plt.plot(xs, ys, label=label)
            if xp is not None and yp is not None:
                plt.plot(xp, yp, 'o', label='data points for ' + label)
        plt.legend()
        plt.show()
        
