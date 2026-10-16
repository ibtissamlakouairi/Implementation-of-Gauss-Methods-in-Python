# Gauss without pivoting
def gauss(A, b):
    n = len(A)
    for k in range(n):
        for i in range(k + 1, n):
            m = A[i][k] / A[k][k]
            for j in range(k, n):
                A[i][j] = A[i][j] - m * A[k][j]
            b[i] = b[i] - m * b[k]
    x = [0.0] * n

    for i in range(n - 1, -1, -1):
        s = 0
        for j in range(i + 1, n):
            s = s + A[i][j] * x[j]
        x[i] = (b[i] - s) / A[i][i]
    return x


# Gauss with partial pivoting
def gauss_partial(A, b):
    n = len(A)
    for k in range(n):
        p = k
        for i in range(k + 1, n):
            if abs(A[i][k]) > abs(A[p][k]):
                p = i
        A[k], A[p] = A[p], A[k]
        b[k], b[p] = b[p], b[k]
        for i in range(k + 1, n):
            m = A[i][k] / A[k][k]
            for j in range(k, n):
                A[i][j] = A[i][j] - m * A[k][j]
            b[i] = b[i] - m * b[k]
    x = [0.0] * n

    for i in range(n - 1, -1, -1):
        s = 0
        for j in range(i + 1, n):
            s = s + A[i][j] * x[j]
        x[i] = (b[i]- s) / A[i][i]
    return x

# Gauss with total pivoting
def gauss_total(A, b):
    n = len(A)
    order = [0, 1, 2]
    for k in range(n):
        p = k
        q = k

        for i in range(k, n):
            for j in range(k, n):
                if abs(A[i][j]) > abs(A[p][q]):
                    p = i
                    q = j

        A[k],A[p]= A[p],A[k]
        b[k],b[p]= b[p],b[k]

        for i in range(n):
            A[i][k],A[i][q] = A[i][q],A[i][k]
        order[k],order[q] = order[q],order[k]

        for i in range(k + 1, n):
            m = A[i][k] / A[k][k]
            for j in range(k, n):
                A[i][j] = A[i][j] - m * A[k][j]
            b[i]= b[i] - m * b[k]
    z = [0.0] * n

    for i in range(n - 1, -1, -1):
        s = 0
        for j in range(i + 1, n):
            s = s + A[i][j] * z[j]
        z[i] = (b[i] - s) / A[i][i]
    x = [0.0] * n

    for i in range(n):
        x[order[i]] = z[i]
    return x

A = [[2, 1, 1], [4, 3, 3], [2, 4, 5]]
b = [8, 20, 21]

print("Without pivoting:")
print(gauss([row[:] for row in A], b[:]))
print("Partial pivoting:")
print(gauss_partial([row[:] for row in A], b[:]))
print("Total pivoting:")
print(gauss_total([row[:] for row in A], b[:]))