import matplotlib.pyplot as plt #install with "pip install matplotlib"

class Stuetzstelle:
    def __init__(self, x, f):
        self.x = x
        self.f = f

    def __str__(self):
        return f'f({self.x}) = {self.f}'

class StuetzstellenListe:
    def __init__(self):
        self.liste = []

    def add(self, sz: Stuetzstelle):
        length = len(self.liste)
        self.liste.append(sz)

    def text(self, i):
        if i >= len(self.liste) or i < 0:
            print('Stuetzstelle nicht vorhanden!')
        return str(self.liste[i])

    def get(self, k, j=0):
        if (k-j) < 0:
            print('k-j darf nicht kleiner als 0 sein!')
            return None
        elif (k-j) >= len(self.liste):
            print('Stuetzstelle nicht vorhanden!')
            return None
        else:
            return self.liste[k-j]


class DreiecksSchema:
    def __init__(self, szs: StuetzstellenListe):
        self.szs = szs
        # init 2d array
        self._k_js = [[None for i in range(j+1)] for j in range(len(szs.liste))]
        n = len(self.szs.liste)
        print(f'Interpolationsproblem mit\nn = {n} \t Stuetzstellen\n')
        if n < 2:
            print('Es müssen mindestens 2 Stützstellen vorhanden sein!')
            return None

    def notation(self, k, j):
        return f'D_{k}_{j}'

    def n_out_of_bounds(self, n):
        return "n muss zwischen 0 und der Anzahl der Stützstellen - 1 liegen!"

    def calc_kj(self, k, j, x=None):
        if j == 0:
            self._k_js[k][j] = self.szs.get(k).f
        elif self._k_js[k][j] is None:
            _k_minus_1_j_minus_1 = self.calc_kj(k-1, j-1, x)
            _k_j_minus_1 = self.calc_kj(k, j-1, x)
            x_k = self.szs.get(k).x
            x_k_minus_j = self.szs.get(k-j).x
            _k_j = self.rekursionsVorschrift(_k_j_minus_1, _k_minus_1_j_minus_1, x_k, x_k_minus_j, x)
            self._k_js[k][j] = _k_j
        return self._k_js[k][j]


    def visualize(self):
        n = len(self.szs.liste) - 1
        # print table
        print('k\\j', end='\t\t')
        for j in range(n+1):
            # print header with j's , no line break
            print(f'j={j}', end='\t\t') # 0, 1, ..., n
        print('\n')
        for k in range(n+1): # 0, 1, ..., n
            # print k's	
            print(f'k={k}, {self.szs.text(k)}', end='\t')
            for j in range(k+1): # 0, 1, ..., k
                print(f'{self.notation(k, j)}={self._k_js[k][j]}', end='\t')
            print('')
        print(f'\n{self.notation(n, n)} = {self._k_js[n][n]}\n')
        print('----------------------------------------\n\n')

class NevilleAitken(DreiecksSchema):
    def __init__(self, szs: StuetzstellenListe):
        super().__init__(szs)

    def notation(self, k, j):
        return f'p_{k}_{j}'

    def rekursionsVorschrift(self, p_k_j_minus_1, p_k_minus_1_j_minus_1, x_k, x_k_minus_j, x):
        p_k_j = p_k_j_minus_1 + (x - x_k) / (x_k - x_k_minus_j) * (p_k_j_minus_1 - p_k_minus_1_j_minus_1)
        return p_k_j

    def interpolate(self, x, n=False, visualize=False):
        if not n:
            n = len(self.szs.liste) - 1
        assert n < len(self.szs.liste) and n >= 0, self.n_out_of_bounds(n)
        p_k_j = self.calc_kj(n, n, x)
        if visualize:
            print(f'Using Neville-Aitken to calculate p_{n}({x})')
            self.visualize()
        return p_k_j

class DividierteDifferezen(DreiecksSchema):
    def __init__(self, szs: StuetzstellenListe):
        super().__init__(szs)

    def notation(self, k, j):
        if j == 0:
            return f'f[x{k}]'
        if j == 1:
            return f'f[x{k}, x{k-1}]'
        return f'f[x{k},...,x{k-j}]'

    def rekursionsVorschrift(self, dd_k_j_minus_1, dd_k_minus_1_j_minus_1, x_k, x_k_minus_j, x):
        # unabhängig von x, x wird nicht verwendet
        return (dd_k_j_minus_1 - dd_k_minus_1_j_minus_1) / (x_k - x_k_minus_j)

    def dds_nth_grade(self, n, visualize=False):
        assert n < len(self.szs.liste) and n >= 0, self.n_out_of_bounds(n)
        dd_k_j = self.calc_kj(n, n)
        if visualize:
            print(f'Calculating Dividierte Differenzen: {self.notation(n, n)}')
            self.visualize()
        return dd_k_j

class Newton:
    def __init__(self, szs: StuetzstellenListe):
        self.szs = szs
        self.dd_calc = DividierteDifferezen(szs)
        self.dds = []

    def basisPolynom(self, n, x, visualize=False):
        product = 1
        if (visualize):
            print(f'Basis polynom w_{0}({x}): {product}')            
        yield product
        product_str = ''
        for i in range(1, n+1): # 1, 2, ..., n
            product *= (x - self.szs.get(i-1).x)
            if (visualize):
                print(f'Basis polynom w_{i}({x}): {product_str} = {product}')            
            yield product

    def interpolate(self, x, n=False, visualize=False):
        if not n:
            n = len(self.szs.liste) - 1
        assert n < len(self.szs.liste) and n >= 0, self.n_out_of_bounds(n)
        local_dds = [] # first grade
        for i in range(0, n):  # 0, 1, ..., n-1, next line n-th grade
            local_dds.append(self.dd_calc.dds_nth_grade(i))  # "caches" the dds       
        local_dds.append(self.dd_calc.dds_nth_grade(n, visualize=visualize))
        self.dds = local_dds
        # using generator function for the basis polynom
        # and using saved list of dds

        # list comprehension
        p_x = sum([dd * product for dd, product in zip(self.dds, self.basisPolynom(n, x))])
        print (f'Newton Darstellung p_{n}({x}) = {p_x}\n')


    def calculate_P_over_interval(self, xstart, xend, step_size):
        x_values = [x for x in range(xstart, xend, step_size)]
        y_values = []
        for x in x_values:
            y_values.append(sum([dd * product for dd, product in zip(self.dds, self.basisPolynom(n, x))]))
       
        # plot
        plt.plot(x_values, y_values)
       



if __name__ == '__main__':
    
    problem = StuetzstellenListe()
    problem.add(Stuetzstelle(0, 2))
    problem.add(Stuetzstelle(1, 1))
    problem.add(Stuetzstelle(3, 2))
    problem.add(Stuetzstelle(4, 4))
    
    aitken = NevilleAitken(problem)
    aitken.interpolate(x=2, n=3, visualize=True)
    

    newton = Newton(problem) # Newton uses Dividierte Differenzen Class
    newton.interpolate(x=2, n=3, visualize=True) # visualize Dividierte Differenzen