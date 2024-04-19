# Kuramoto_sivashinsky

**This repo will cheese the KS equation. **  

KS_Fourier_NR V2 == uses newton raphson to find periodic orbit listed in this paper: 
```diff
+ Spatiotemporal chaos in terms of unstable recurrent patterns
- F Christiansen1, P Cvitanovic1 and V Putkaradze1
- Published under licence by IOP Publishing Ltd   
```   

KS_Fourier_NR V2 TW == Uses newton raphson to find Travelling wave listed in this paper
```diff
+ On the state space geometry of the Kuramoto-Sivashinsky flow in a periodic domain
- Predrag Cvitanović, Ruslan L. Davidchack, Evangelos Siminos
```
Finite difference solver:
Finite_difference_KS_Newton_raphson: Use newton raphson to find equilibrium using finite difference discretization.

Spectral method solvers
Fourier_KS_Newton_Raphson: Use Newton raphson to find equilibrium of KS equation
KS_Fourier_NR V2: Use newton raphson to find periodic orbit of KS equation. 
KS_Fourier_NR V2 TW: Use newton raphson to find Travelling wave of KS equation

Dedalus solver:
KS_GD_Using_Adjoint_Operator:  Folder contains the use of adjoint operator to find equilibriums of KS equation
The method is called ajoint descent which is similar to gradient descent, but here we have made a second PDE with 
adjoint operator which in the fictitious time evolution, reduce the functional norm, gradually taking us to an equilibirium.
This method works faster than regular gradient descent method.
