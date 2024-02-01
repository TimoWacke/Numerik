'''
    Stuetzstelle Class, Stuetzstellenlist Class and TriangularRecursive
Schema Class are useful for
    the recursive scheme...

    With the DividedDifferences
 Class that is a subclass of TriangularRecursive
Schema
    and the Newton Class we have the tools to interpolate a given function with not only
    Newton Interpolation but also with minimal effort with Hermite Interpolation or Tschebyscheff Interpolation

    These classes have been used in a previous exercise already
'''


class SupportingPoint:
    def __init__(self, x, f):
        self.x = x
        self.f = f

    def __str__(self):
        return f'f({self.x}) = {self.f:.2f}'

class SupportPointList:
    def __init__(self):
        self.list = []

    def add(self, sp: SupportingPoint):
        self.list.append(sp)

    def get(self, k, j=0):
        if (k-j) < 0:
            print('k-j should not be negative!')
            return None
        elif (k-j) >= len(self.list):
            print('k-j should not be greater than the length of the list!')
            return None
        else:
            return self.list[k-j]
        
    def get_xlist(self):
        return [sp.x for sp in self.list]
    
    def get_ylist(self):
        return [sp.f for sp in self.list]

class TriangularRecursiveSchema:
    def __init__(self, spl: SupportPointList):
        self.spl = spl
        if len(self.spl.list) > 0:
            # init 2d array
            self._k_js = [[None for i in range(j+1)] for j in range(len(spl.list))]
        else:
            self._k_js = []
        
    def expand_size_to(self, n):
        # add empty rows to the 2d array
        for i in range(len(self._k_js), n):
            self._k_js.append([None for i in range(i+1)])

    def notation(self, k, j):
        return f'D_{k}_{j}'

    def calc_kj(self, k, j):
        if j == 0:
            self._k_js[k][j] = self.spl.get(k).f
        elif self._k_js[k][j] is None:
            _k_minus_1_j_minus_1 = self.calc_kj(k-1, j-1)
            _k_j_minus_1 = self.calc_kj(k, j-1)
            x_k = self.spl.get(k).x
            x_k_minus_j = self.spl.get(k-j).x
            _k_j = self.recursion_rule(_k_j_minus_1, _k_minus_1_j_minus_1, x_k, x_k_minus_j)
            assert _k_j is not None
            self._k_js[k][j] = _k_j
        return self._k_js[k][j]


    def visualize(self):
        n = len(self.spl.list) - 1
        # print table
        print('k\\j', end='\t\t')
        for j in range(n+1):
            # print header with j's , no line break
            print(f'j={j}', end='\t\t') # 0, 1, ..., n
        print('\n')
        for k in range(n+1): # 0, 1, ..., n
            # print k's	
            print(f'k={k}, {self.spl.get(k).x:.2f}', end='\t')
            for j in range(k+1): # 0, 1, ..., k
                print(f'{self.notation(k, j)}={self._k_js[k][j]:.2f}', end='\t')
            print('')
        print(f'\n{self.notation(n, n)} = {self._k_js[n][n]}\n')
        print('----------------------------------------\n\n')


class DividedDifferences(TriangularRecursiveSchema):
    def __init__(self, spl):
        super().__init__(spl)

    def notation(self, k, j):
        if j == 0:
            return f'f[x{k}]'
        if j == 1:
            return f'f[x{k}, x{k-1}]'
        return f'f[x{k}..x{k-j}]'

    def recursion_rule(self, dd_k_j_minus_1, dd_k_minus_1_j_minus_1, x_k, x_k_minus_j):
        dd_k_j = (dd_k_j_minus_1 - dd_k_minus_1_j_minus_1) / (x_k - x_k_minus_j)
        return dd_k_j

    def dds_nth_grade(self, n, visualize=False):
        assert n < len(self.spl.list) and n >= 0, self.n_out_of_bounds(n)
        dd_k_j = self.calc_kj(n, n)
        if visualize:
            print(f'Calculating divided difference: {self.notation(n, n)}')
            self.visualize()
        return dd_k_j
   