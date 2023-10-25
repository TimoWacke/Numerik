

'''
    Definition einer Rekursiven Funktion um die Dividierten Differenz f[x_(k-j), ..., x_k] zu berechnen
    Es werden die x_i und f_i als Listen übergeben, sowie k und j als Indizes
'''
def dividierte_differenzen(xs, fs, k, j):
    assert k >= 0 and k < len(xs), 'k muss zwischen 0 und der Anzahl der Stützstellen - 1 liegen!'
    assert j >= 0 and j <= k, 'j muss zwischen 0 und k liegen!'
    assert len(xs) == len(fs), 'xs und fs müssen gleich lang sein!'

    # Abbruchbedingung
    if j == 0:
        dd_k_j = fs[k]
    else:
        dd_k_j_minus_1 = dividierte_differenzen(xs, fs, k, j-1) # rekursiver Aufruf von sich selbst
        dd_k_minus_1_j_minus_1 = dividierte_differenzen(xs, fs, k-1, j-1) # rekursiver Aufruf von sich selbst
        x_k = xs[k]
        x_k_minus_j = xs[k-j]
        dd_k_j = (dd_k_j_minus_1 - dd_k_minus_1_j_minus_1) / (x_k - x_k_minus_j)
    return dd_k_j


# Aufruf mit den Stützstellen aus Aufgabe 2

stuetzstellen_x = [0, 1, 3, 4]
stuetzstellen_f = [2, 1, 2, 4]
print("Dividierte Differenz f[x_0, ..., x_3] für die Stützstellen aus Aufgabe 2:")
print(dividierte_differenzen(stuetzstellen_x, stuetzstellen_f, 3, 3),)