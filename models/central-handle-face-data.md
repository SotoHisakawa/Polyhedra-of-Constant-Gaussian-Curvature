# Face data for the central handle

The development in Appendix B.2 shows one face of each nonrectangular type. The construction uses a regular prism of base side length $2$ and half-height $h$. The paper treats $4\le g\le7$ with $h>0$ sufficiently small. The genus-five model used for the development has $h=1/3$. The function `build_central_handle` generates its coordinates and faces. [Model coordinates and faces](embedded_orientable_g5_v28_central_handle.json) · [Python construction](../scripts/embedded_models.py)

Write $n=g-1$, $a=\pi/n$, $b=a+\pi/14$, $\kappa=2\pi/7$, $\rho=\cot a$, $\tau=\cot b$, and $q=\rho-\tau$. The window displacement $u$, window half-width $y_0$ at the waist, and inward displacement $x_0$ of the middle ring are

$$
u=\frac{h\sin b\sin\kappa}{\sin b+\cos\kappa},\qquad
y_0=1-u,\qquad
x_0=\frac{h\sin(\pi/7)}{\sqrt{\sin^2a-\sin^2(\pi/7)}}.
$$

The construction requires $0<u<1$ and $0<x_0<q$. The lengths of the sloping edges of the window and central tunnel are

$$
L_w=\sqrt{u^2\csc^2 b+h^2},\qquad L_i=\sqrt{x_0^2+h^2}.
$$

At $v_1$, the pentagon has angle $\psi=\arctan(h\sin b/u)$ and the window trapezoid has angle $\zeta=\arccos(\sin b\cos\psi)$. They satisfy $\psi+\zeta=\pi-\kappa=5\pi/7$.

| Face | Vertices in cyclic order | Successive edge lengths | Interior angles in the same order |
|---|---|---|---|
| Window pentagon | $v_1,v_3,v_4,v_2,v_5$ | $\csc b,\ 2h,\ \csc b,\ L_w,\ L_w$ | $\psi,\ \pi/2,\ \pi/2,\ \psi,\ 2\pi-2\psi$ |
| Window trapezoid | $v_1,v_5,v_6,v_{1,1}$ | $L_w,\ 2y_0,\ L_w,\ 2$ | $\zeta,\ \pi-\zeta,\ \pi-\zeta,\ \zeta$ |
| Top triangle | $v_1,v_3,v_{3,-1}$ | $\csc b,\ 2q\sin a,\ \csc b$ | $\pi/7,\ 3\pi/7,\ 3\pi/7$ |
| Central tunnel trapezoid | $v_{3,-1},v_{7,-1},v_7,v_3$ | $L_i,\ 2(q-x_0)\sin a,\ L_i,\ 2q\sin a$ | $5\pi/14,\ 9\pi/14,\ 9\pi/14,\ 5\pi/14$ |

In each row the first length joins the first and second vertices. Angles are in radians. The angle at $v_5$ in the pentagon is reflex.

### Indices in the genus-five JSON file

The vertex $v_{i,k}$ is $v_i$ rotated through $2\pi k/n$ about the vertical axis, with $v_i=v_{i,0}$. Its zero-based index in the `vertices` array is $7(k\bmod n)+i-1$. The following indices follow the cyclic orders used above. A cycle in the `faces` array may be reversed.

| Face | Index in `faces` | Indices in `vertices` |
|---|---:|---|
| Window pentagon | 5 | `0, 2, 3, 1, 4` |
| Window trapezoid | 3 | `0, 4, 5, 7` |
| Top triangle | 1 | `0, 2, 23` |
| Central tunnel trapezoid | 34 | `23, 27, 6, 2` |
