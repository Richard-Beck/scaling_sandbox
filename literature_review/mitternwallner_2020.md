Mitterwallner et al., 2020, “Non-Markovian data-driven modeling of single-cell motility,” Physical Review E 101, 032408. DOI: 10.1103/PhysRevE.101.032408

Mitterwallner et al. model the 1D tangential cell velocity with a generalized Langevin equation

$$
\dot v(t)
=
-\int_{-\infty}^{t}\Gamma_v(t-s)\,v(s)\,ds
+F_R(t),
$$

ignoring an external positional potential for the circular-track case. The random force can itself be colored, with

$$
\langle F_R(t)F_R(0)\rangle=\Gamma_R(t).
$$

Because the trajectory alone cannot uniquely identify both \(\Gamma_v\) and \(\Gamma_R\), they introduce a **substitute kernel**

$$
\Gamma(t)\equiv \Gamma_v(t)=\frac{\Gamma_R(t)}{B},
\qquad
B=\langle v^2\rangle .
$$

([APS Journals][1])

What the data give them is the particularly simple kernel

$$
\boxed{
\Gamma(t)=2a\,\delta(t)+b\,e^{-t/\tau_m},
\qquad b<0 .
}
$$

So the velocity feels ordinary instantaneous damping plus a **delayed negative-friction term**: previous motion pushes the present velocity in the same direction, with memory dying away over \(\tau_m\). A standard OU/persistent-random-walk model is basically the \(b=0\) limit. ([APS Journals][1])

For their fit to the **cell-averaged MDA-MB-436 data**, they report

$$
\boxed{
a=0.15~{\rm min}^{-1},\qquad
b=-0.023~{\rm min}^{-2},\qquad
\tau_m=6.02~{\rm min}
}
$$

and

$$
\boxed{
B=\langle v^2\rangle=0.14~\mu{\rm m}^2/{\rm min}^2
}
$$

so an RMS velocity scale is

$$
\sqrt B \approx 0.374~\mu{\rm m}/{\rm min}.
$$

Those parameters generate two emergent relaxation times,

$$
\tau_b=3.23~{\rm min},
\qquad
\boxed{\tau_p=164~{\rm min}},
$$

and an effective long-time diffusivity

$$
D=12.7~\mu{\rm m}^2/{\rm min}.
$$

([APS Journals][1])

That huge separation is the interesting bit: the **actual kernel only remembers velocity for \(\sim6\) min**, yet the near-cancellation

$$
a \simeq |b|\tau_m
$$

between positive instantaneous friction and integrated negative delayed friction creates an emergent persistence time of roughly **2.7 hours**. Across individual cells they actually find this relation \(a\approx |b|\tau_m\) quite strongly. ([APS Journals][1])

So, written directly as a simulation equation using their average parameters:

$$
\boxed{
\dot v(t)
=
-0.15\,v(t)
+
0.023\int_{-\infty}^{t}
e^{-(t-s)/6.02}\,v(s)\,ds
+
F_R(t)
}
$$

with consistent minute/\(\mu\)m units.

And computationally you don't need to retain the entire history: because the kernel is exponential, introduce

$$
h(t)=\int_{-\infty}^{t}e^{-(t-s)/\tau_m}v(s)\,ds,
$$

so

$$
\dot h=v-\frac{h}{\tau_m},
$$

giving a **2-state Markovian system**

$$
\dot v=-av-bh+F_R(t),
\qquad
\dot h=v-\frac{h}{\tau_m},
$$

where remember \(b<0\), so \(-bh\) is positive feedback. That's probably the most useful form for what you're building.

[1]: https://journals.aps.org/pre/pdf/10.1103/PhysRevE.101.032408 "Non-Markovian data-driven modeling of single-cell motility"
