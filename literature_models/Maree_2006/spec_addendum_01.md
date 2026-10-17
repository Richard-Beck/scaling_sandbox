Yes. **The issue is real, and I would amend my previous specification.** The literal combination “\(a,A\) in nm² + \(\lambda=10^{-9}\)” cannot be the implementation that generated the published simulations. Your diagnostic demonstrates that numerically, and the paper is internally dimensionally inconsistent on exactly this point.

The strongest evidence is:

* Marée defines

  $$
  H=H_{\rm contact}+\lambda(a-A)^2
  $$

  and describes \(a,A\) as cell areas. ([ResearchGate][1])
* Table 3 simultaneously gives \(T=0.008\,{\rm nm^{-1}}\), \(H_b=0.046\,{\rm nm^{-1}}\), \(J=7.5\times10^{-4}\,{\rm nm^{-1}}\), \(\lambda=10^{-9}\,{\rm nm^{-3}}\), and \(A=3\times10^8\,{\rm nm^2}\). ([ResearchGate][1])
* But

  $$
  [\lambda(a-A)^2]
  ={\rm nm^{-3}}\times({\rm nm^2})^2
  ={\rm nm},
  $$

  whereas it has to be commensurate with \(T,H_b,J,P\), i.e. nm\(^{-1}\). So the published units cannot all be correct.
* The authors explicitly say the same target area is

  $$
  3\times10^8\,{\rm nm^2}\simeq3\times10^4\ {\rm lattice\ sites}.
  $$

  ([ResearchGate][1])
* They also say \(\lambda\) has only a **“limited effect on the dynamics”** and is chosen merely to keep cell-size fluctuations within a few percent. Your literal implementation instead makes \(\lambda\) overwhelmingly dominate the dynamics and freezes the membrane. That directly contradicts their description. ([ResearchGate][1])
* Their standard full model demonstrably does undergo large numbers of boundary updates and moves persistently at about \(90\,{\rm nm/s}\). ([ResearchGate][1])

There is an additional piece of evidence I think is quite important. I checked the **Tissue Simulation Toolkit code lineage from the same Utrecht CPM ecosystem**. The current source says explicitly that it derives from a CPM implementation written around 1995, with copyright extending through 2006. Its cell `Area()` and `TargetArea()` are integer quantities, and its area-energy difference for an extension is implemented as

$$
\Delta H_A
=
\lambda\left[1+2(A_{\rm current}-A_{\rm target})\right],
$$

i.e. **one accepted pixel changes area by exactly 1**, not by \(10^4\ {\rm nm^2}\).

That isn't proof that the keratocyte program used precisely this source file, but given the date, group and CPM lineage, it's substantially stronger evidence than merely saying “CPMs usually count pixels.”

### What I think happened

I now think the paper **mixed two descriptions of area**:

$$
A_{\rm CPM}\simeq30,000\ {\rm sites}
$$

for the actual Hamiltonian computation, and

$$
A_{\rm physical}=A_{\rm CPM}\Delta x^2
\simeq 3\times10^8\ {\rm nm^2}
$$

for reporting the biological size.

Somewhere in converting the CPM parameters to physical units for the paper, the units attached to \(\lambda\) became inconsistent. There may also have been an incorrect dimensional conversion in the paper itself. There is already another suspicious unit in Table 3: the printed units for \(\xi\) don't match Eq. 24 dimensionally either.

However, there is one complication: **simply using \(\lambda=10^{-9}\) with site-count area may also not be the final answer.** At target area it gives

$$
\Delta H_A = 10^{-9},
$$

which is utterly negligible compared with the \(\sim10^{-3}\) contact contribution. Even hundreds of pixels from target it remains tiny. So although the lattice-site interpretation fixes your catastrophic \(0.1\) barrier, we still have to establish whether the numerical \(\lambda=10^{-9}\) was itself rescaled before entering their CPM code.

That is now the key unresolved reconstruction issue.

### Why I would **not** make \(10^{-10}\) your preferred fix

Your proposed \(10^{-10}\) physical-area coefficient is an excellent **diagnostic**, but it's not a principled conversion.

If

$$
a_{\rm phys}=N\,A_{\rm pix},
\qquad A_{\rm pix}=10^4\,{\rm nm^2},
$$

then

$$
\lambda_{\rm phys}(a_{\rm phys}-A_{\rm phys})^2
=
\lambda_{\rm lattice}(N-N_0)^2
$$

requires

$$
\boxed{
\lambda_{\rm phys}
=
\frac{\lambda_{\rm lattice}}{A_{\rm pix}^2}
}
$$

and hence a factor

$$
A_{\rm pix}^2=10^8.
$$

So if the literal numerical CPM coefficient really were \(10^{-9}\) acting on site count, its physical-area equivalent would be

$$
\lambda_{\rm phys}=10^{-17}\ {\rm nm^{-5}},
$$

**not \(10^{-10}\)**.

Conversely, using \(10^{-10}\) with nm² areas corresponds to an enormous lattice-area coefficient:

$$
\lambda_{\rm lattice}=10^{-10}(10^4)^2=10^{-2}.
$$

That can certainly make the simulation move, but it describes a very stiff area constraint quite different from “\(\lambda=10^{-9}\)” in lattice units.

### What I would implement now

I would retain three clearly separated variants rather than silently retuning anything:

1. **`published_literal`** — your present control:

   $$
   a=N10^4\,{\rm nm^2},\quad
   A=3\times10^8\,{\rm nm^2},\quad
   \lambda=10^{-9}.
   $$

   It freezes. Keep this result because it proves the paper cannot be interpreted literally.

2. **`cpm_lattice_area` — my new primary historical reconstruction hypothesis:**

   $$
   a=N,\qquad A=30,000,
   $$

   with an initially separate parameter \(\lambda_{\rm CPM}\). Start by testing the reported numerical value \(10^{-9}\), but don't assume yet that it is the final code value.

   The exact incremental term should then be

   $$
   \Delta H_A^{\rm ext}
   =
   \lambda_{\rm CPM}[2(N-N_0)+1],
   $$

   $$
   \Delta H_A^{\rm ret}
   =
   \lambda_{\rm CPM}[-2(N-N_0)+1].
   $$

3. **`physical_area_equivalent`**, only as a unit-converted implementation of #2:

   $$
   \lambda_{\rm phys}
   =
   \frac{\lambda_{\rm CPM}}{10^8}\ {\rm nm^{-5}}.
   $$

   This should produce **identical copy energies** to #2.

That #2/#3 equivalence itself should become an automated unit test.

Then determine \(\lambda_{\rm CPM}\) using the paper's actual stated constraint: **cell-size fluctuations remain within a few percent while \(\lambda\) has limited influence on motility**. ([ResearchGate][1]) I would sweep \(\lambda_{\rm CPM}\) logarithmically and measure:

$$
\mathrm{CV}(A),\qquad
\frac{\langle A\rangle-A_0}{A_0},
\qquad
v_{\rm COM},
\qquad
\text{copy rate}.
$$

The successful range should retain plentiful membrane copying, maintain area around \(30{,}000\) sites, and then still permit the published \(90\,{\rm nm/s}\)-ish motile phenotype after polarization. That's much more informative than asking only where copying first becomes possible.

One other useful constraint is their force–velocity calibration. Their derivation of \(T=0.008\) and \(H_b=0.046\) treats protrusion/retraction as readily occurring stochastic events governed by pushing barbed ends and membrane yield. ([ResearchGate][1]) A hidden \(+0.1\) area penalty on every extension would completely invalidate that calibration. So your result doesn't merely look suspicious: **it is incompatible with the derivation used to set two of the model's central mechanical parameters.**

I therefore think your diagnostic has found a genuine error in my reproduction spec. The important revision is **not “change \(\lambda\) to \(10^{-10}\)”**, but **“do not feed physical nm² directly into Eq. 22 with the printed \(\lambda\); reconstruct the CPM's internal area convention first.”**

I also tried to inspect your `BOUNDARY_COPY_DIAGNOSTIC.md` through the connected `scaling_sandbox` repo, but that path isn't currently visible on the connected/default branch, so the assessment of those run numbers is based on what you pasted here. The paper/TST evidence above is independently accessible.

[1]: https://www.researchgate.net/publication/6989407_Polarization_and_Movement_of_Keratocytes_A_Multiscale_Modelling_Approach "(PDF) Polarization and Movement of Keratocytes: A Multiscale Modelling Approach"
