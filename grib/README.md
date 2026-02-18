# OPCONTI Grib

📜 _Redacted notes, reports, papers and presentations ..._

## Details

- [01_notes](./01_notes) : Liudi's notes from the discussion in Linz & Oberwolfach

> 📣 PDF files are tracked by `git`, so one doesn't need to compile the `latex` sources to view its content. But in case of merge conflicts, always **solve first the conflict on the `latex` file(s), then re-compile the related PDF file(s), and commit afterward**. 
> Never try to solve a conflict on a binary PDF file ... 😇

## Miscellaneous

_Short summary on the late discussions ..._

> Meeting 18/02/2026 (GC, LL, TL, TV)
  see ![notes from Gabriele](PhotoGabriele.png)

> Meeting 20/01/2026 (LL, GC, TL)

- idea of using a coarse time grid when decomposing the variables
- first tests from LL on the domain decomposition
- objectives : write down the decomposition on the continuous level (theory behind the code)
  - stick to linear quadratic
  - solve globally a using a "local" right-hand-side 

> [GC]

A) I would like to consider at least three test cases:

- a linear quadratic one (as written in the lecture notes)
- a bilinear one (Tom, one like the problem with Dominique, without stochasticity of course)
- a fully nonlinear one (we can think about that)

The point is that, a linear quadratic problem will lead to a quadratic minimization for which we know that there is an acceleration (see the [proceeding with Tom](https://arxiv.org/pdf/2409.01032)). Indeed, the effect of all usual optimal control parameters need to be investigated. Moreover, the linear quadratic one could be the first "easy" test to analyze the behavior.

B) In our preliminary tests, I would consider a Crank-Nicolson type discretization for the forward problem and a discretize-before-optimize approach for the optimality system.
This implicit method rules out numerical instability and guarantees norm preserving of the test case 2).

> [TL]

A) for the third test case I can suggest the Shallow Water Equation problem [as used here](https://arxiv.org/pdf/2404.05556) : already have the code developed by the PhD student (is in Python tho 😉).

B) for the time-integration scheme, Crank-Nicolson seems fine although Backward Euler would provide the same numerical stability if that's what we want. 
But if we use an implicit time scheme, we still have to build a Newton based solver for the ODE RHS at some point. E.g if we consider :

$$
\frac{\partial {\bf y}}{\partial t} = f({\bf y}, {\bf u})
$$

with ${\bf y}$ the state and ${\bf u}$ the control.
We need a routine computing the solution ${\bf y}$ of 

$$
{\bf y} - \alpha f({\bf y}, {\bf u}) = {\bf r}
$$

for whatever $\alpha$ and ${\bf r}$ we plug in. 
Usually, there are some generic non-linear solver in Matlab/Python 
based on MINPACK allowing to solve that using only evaluations of $f$ (`fsolve`).
I would start by using those to implement the first test cases 
$\Rightarrow$ just need the definition of $f({\bf y}, {\bf u})$
