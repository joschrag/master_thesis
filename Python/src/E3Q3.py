import numpy as np
import numpy.linalg as la
import sympy as sp
from sympy import Eq, simplify
from sympy.matrices import Matrix

x, y, z = sp.symbols("x, y, z")

if __name__ == "__main__":
    m = np.array(
        [
            [1, 2, 1, 1, 1, 0, 1, 1, 1, -10],
            [1, 1, 3, 0, 1, -1, 7, 0, 0, -10],
            [1, 1, 1, 0, 0, 1, -15, 0, 0, -10],
        ],
        dtype=float,
    )

    A = -np.transpose(np.array([m[:, 1], m[:, 2], m[:, 5]]))
    print(A, la.matrix_rank(A), m[:, 1], sep="\n")
    if la.matrix_rank(A) == 3:
        Q = Matrix(la.inv(A))
        P = Matrix(
            [
                [
                    m[0, 3] * x + m[0, 7],
                    m[0, 4] * x + m[0, 8],
                    m[0, 0] * x**2 + m[0, 6] * x + m[0, 9],
                ],
                [
                    m[1, 3] * x + m[1, 7],
                    m[1, 4] * x + m[1, 8],
                    m[1, 0] * x**2 + m[1, 6] * x + m[1, 9],
                ],
                [
                    m[2, 3] * x + m[2, 7],
                    m[2, 4] * x + m[2, 8],
                    m[2, 0] * x**2 + m[2, 6] * x + m[2, 9],
                ],
            ],
        )
        sp.pprint(P)
        P2 = Q * P
        sp.pprint(P2)
        lin_vars = Matrix([y, z, 1])
        quad_1 = (P2.row(0) * lin_vars)[0]
        quad_2 = (P2.row(1) * lin_vars)[0]
        mixed = (P2.row(2) * lin_vars)[0]

        print(quad_1, quad_2, mixed, sep="\n")
        identities = np.array(
            [
                Eq(
                    (quad_1 * lin_vars[1]),
                    (mixed * lin_vars[0]),
                    evaluation=False,
                ),
                Eq(
                    (mixed * lin_vars[1]),
                    (quad_2 * lin_vars[0]),
                    evaluation=False,
                ),
                Eq((mixed * mixed), (quad_1 * quad_2), evaluation=False),
            ]
        )
        substitutions = {
            lin_vars[0] ** 2: quad_1,
            lin_vars[1] ** 2: quad_2,
            lin_vars[0] * lin_vars[1]: mixed,
        }
        for i in range(2):
            identities = np.array(
                [
                    simplify(identity.expand(basic=True).subs(substitutions))
                    for identity in identities
                ]
            )
            print(f"{i+1}:{identities=}")

        identities = np.array(
            [
                Eq(
                    (identity.lhs - identity.rhs).expand(basic=True),
                    0,
                )
                for identity in identities
            ]
        )
        print(f"3:{[i.lhs for i in identities]=}")
        identities = np.array([i.lhs for i in identities])
        print("#############################################")
        print(identities)
        print("#############################################")
        A, b = sp.linear_eq_to_matrix(identities, [y, z])
        A = A.col_insert(2, -b)
        sp.pprint(A)
        poly: sp.Poly = sp.Poly(A.det())
        coeffs = np.array(poly.all_coeffs(), dtype=np.float64)[::-1]
        f = np.polynomial.Polynomial(coeffs)
        x_sol = f.roots()
        x_sol = x_sol.real[abs(x_sol.imag) < 1e-5]
        x_0 = x_sol[0]
        for x_0 in x_sol:
            M = np.array(A.subs({x: x_0})).astype(np.float64)

            _, S, V = la.svd(M)
            print(S, V, sep="\n")
            index = np.argmin(S)
            y_0 = V[index, 0] / V[index, 2]
            z_0 = V[index, 1] / V[index, 2]
            print(f"{x_0=}, {y_0=}, {z_0=}")
        print(
            m
            @ np.array(
                [
                    [x_0**2],
                    [y_0**2],
                    [z_0**2],
                    [x_0 * y_0],
                    [x_0 * z_0],
                    [y_0 * z_0],
                    [x_0],
                    [y_0],
                    [z_0],
                    [1],
                ]
            )
        )
