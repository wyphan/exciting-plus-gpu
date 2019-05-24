#!/bin/bash
cat << "__EOF__" > tasks.tex
\begin{longtable}{lX}
   -1 & Write out the version number of the code. \\
   0 & Ground state run starting from the atomic densities. \\
   1 & Resumption of ground state run using density in {\tt STATE.OUT}. \\
   2 & Structural optimisation run starting from the atomic densities, with
    atomic positions written to {\tt GEOMETRY.OUT}. \\
   3 & Resumption of structural optimisation run using density in
    {\tt STATE.OUT} but with positions from {\tt elk.in}. \\
   5 & Ground state Hartree-Fock run. \\
   10 & Total, partial and interstitial density of states (DOS). \\
   14 & Plots the smooth Dirac delta and Heaviside step functions used by the
        code to calculate occupancies. \\
   15 & Output ${\bf L}$, ${\bf S}$ and ${\bf J}$ total expectation values. \\
   16 & Output ${\bf L}$, ${\bf S}$ and ${\bf J}$ expectation values for each
        {\bf k}-point and state in {\tt kstlist}. \\
   20 & Band structure plot. \\
   21 & Band structure plot which includes angular momentum characters for
    every atom. \\
   25 & Compute the effective mass tensor at the {\bf k}-point given by
    {\tt vklem}. \\
   31, 32, 33 & 1/2/3D charge density plot. \\
   41, 42, 43 & 1/2/3D exchange-correlation and Coulomb potential plots. \\
   51, 52, 53 & 1/2/3D electron localisation function (ELF) plot. \\
   61, 62, 63 & 1/2/3D wavefunction plot:
    $\left|\Psi_{i{\bf k}}({\bf r})\right|^2$. \\
   72, 73 & 2/3D plot of magnetisation vector field, ${\bf m}({\bf r})$. \\
   82, 83 & 2/3D plot of exchange-correlation magnetic vector field,
    ${\bf B}_{\rm xc}({\bf r})$. \\
   91, 92, 93 & 1/2/3D plot of $\nabla\cdot{\bf B}_{\rm xc}({\bf r})$. \\
   100 & 3D Fermi surface plot using the scalar product
    $p({\bf k})=\Pi_i(\epsilon_{i{\bf k}}-\epsilon_{\rm F})$. \\
   101 & 3D Fermi surface plot using separate bands (minus the Fermi
    energy). \\
   110 & Calculation of M\"{o}ssbauer contact charge densities and magnetic
    fields at the nuclear sites. \\
   115 & Calculation of the electric field gradient (EFG) at the nuclear
    sites. \\
   120 & Output of the momentum matrix elements
    $\langle\Psi_{i{\bf k}}|-i\nabla|\Psi_{j{\bf k}}\rangle$. \\
   121 & Linear optical response tensor. \\
   122 & Magneto optical Kerr effect (MOKE) angle. \\
   130 & Output matrix elements of the type
    $\langle\Psi_{i{\bf k+q}}|\exp[i({\bf G+q})\cdot{\bf r}]|
    \Psi_{j{\bf k}}\rangle$. \\
   140 & Energy loss near edge structure (ELNES). \\
   142, 143 & 2/3D plot of the electric field
    ${\bf E}({\bf r})\equiv\nabla V_{\rm C}({\bf r})$. \\
   152, 153 & 2/3D plot of
    ${\bf m}({\bf r})\times{\bf B}_{\rm xc}({\bf r})$. \\
   162 & Scanning-tunneling microscopy (STM) image. \\
   170 & (undocumented) Calls {\tt writewfpw}. \\
   172 & (undocumented) Calls {\tt dielectricq}. \\
   175 & (undocumented) Calls {\tt yambo}. \\
   190 & Write the atomic geometry to file for plotting with {\sf XCrySDen}
    and {\sf V\_Sim}. \\
   200 & Calculation of dynamical matrices on a {\bf q}-point set defined by
    {\tt ngridq}. \\
   210 & Phonon density of states. \\
   220 & Phonon dispersion plot. \\
   230 & Phonon frequencies and eigenvectors for an arbitrary
    ${\bf q}$-point. \\
   240 & Generate the ${\bf q}$-dependent phonon linewidths and electron-phonon
    coupling constants and write them to file. \\
   245 & Phonon linewidths plot. \\
   250 & Eliashberg function $\alpha^2F(\omega)$, electron-phonon coupling
    constant $\lambda$, and the McMillan-Allen-Dynes critical temperature
    $T_c$. \\
   300 & Reduced density matrix functional theory (RDMFT) calculation. \\
   400 & Calculation of tensor moments and corresponding LDA+U Hartree-Fock
    energy contributions. \\
   700 & (undocumented) Calls {\tt sic\_gndstate}. \\
   701 & (undocumented) Calls {\tt sic\_main}. \\
   800 & (undocumented) Calls {\tt response}. \\
   801 & (undocumented) Calls {\tt crpa}. See {\tt CRPA-Calculation.pdf} for details. \\
   802 & (undocumented) Calls {\tt gwmain}. \\
   804 & (undocumented) Calls {\tt genscell}. \\
   805 & (undocumented) Calls {\tt genwfdrc}. \\
   806 & (undocumented) Calls {\tt writebz}. \\
   807, 808 & (undocumented) Calls {\tt writewann}. \\
   811 & (undocumented) Calls {\tt dosrlm}. This is preferred over task 10. \\
   822 & (undocumented) Calls {\tt bandrlm}. This is preferred over task 20. \\
   861, 862, 863 & (undocumented) Calls {\tt wann\_plot}, {\tt wann\_plot\_2d}, {\tt wann\_plot\_3d}, respectively. \\
   880 & (undocumented) Calls {\tt test\_xc}. \\
   881 & (undocumented) Calls {\tt test\_bloch\_wf}. \\
   882 & (undocumented) Calls {\tt writewf}.
\end{longtable}
__EOF__
