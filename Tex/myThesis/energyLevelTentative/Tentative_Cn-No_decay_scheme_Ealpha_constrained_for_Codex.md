# Tentative Cn–No decay scheme: $E_\alpha$-constrained qualitative level heights for Codex

> **Purpose:** This is the corrected drawing specification for Codex.  
> Unlike the previous version, the **relative level heights must be constrained by the observed $E_\alpha$ values** wherever the data allow it.  
> All numerical relations below use the non-relativistic recoil correction
>
> $$
> Q_\alpha \approx E_\alpha\frac{A}{A-4},
> $$
>
> and ignore the small electron-screening correction.
>
> **Important:** the scheme is still tentative. Some of the four $^{265}$Sg $\rightarrow$ $^{261}$Rf arrows are alternative endpoint assignments of the same observed parent $\alpha$ line, not four independently measured $\alpha$ energies.

---

# 1. Experimental $E_\alpha$ values used

Use the following observed energies:

| Parent / branch | $E_\alpha$ |
|---|---:|
| $^{277}$Cn working branch | $11.11$ MeV |
| $^{273}$Ds$_a$ | $10.858$ MeV |
| $^{273}$Ds$_c$ | $11.165$ MeV |
| $^{269}$Hs high-energy group | $9.13$ MeV |
| $^{269}$Hs low-energy group | $8.95$ MeV |
| $^{265}$Sg$_a$ | $8.84$ MeV |
| $^{265}$Sg$_b$ | $8.69$ MeV |
| $^{261}$Rf$_a$ | $8.28$ MeV |
| $^{261}$Rf$_b$ | $8.51$ MeV |

The rejected $^{277}$Cn $\rightarrow {}^{273}$Ds $l_\alpha=5$ branch must not be used.

---

# 2. General energy relation

For two decay paths of the same parent/daughter mass pair,

$$
Q_{\alpha,1}-Q_{\alpha,2}
=
\left(E_{p,1}-E_{p,2}\right)
-
\left(E_{d,1}-E_{d,2}\right).
$$

Therefore the difference of measured alpha energies constrains a **difference of parent and daughter excitation-energy splittings**, rather than fixing both splittings independently.

---

# 3. $^{277}$Cn $\rightarrow {}^{273}$Ds

Use only

$$
^{277}{\rm Cn}(3/2^+[611])
\xrightarrow[\ 11.11\ {\rm MeV}\ ]{\alpha}
^{273}{\rm Ds}(3/2^+[611]).
$$

Do not draw a Cn branch into $^{273}$Ds$(13/2^-[716])$.

Because only one Cn $\rightarrow$ Ds working branch is retained, the Cn alpha energy by itself gives **no relative-energy constraint** between the two Ds states.

---

# 4. $^{273}$Ds $\rightarrow {}^{269}$Hs

Assignments:

$$
^{273}{\rm Ds}_c(3/2^+[611])
\xrightarrow[\ 11.165\ {\rm MeV}\ ]{\alpha}
^{269}{\rm Hs}
\left[
3/2^+\left(1/2^+[620]\right)
\right],
$$

and

$$
^{273}{\rm Ds}_a(13/2^-[716])
\xrightarrow[\ 10.858\ {\rm MeV}\ ]{\alpha}
^{269}{\rm Hs}(11/2^-[725]).
$$

The recoil-corrected alpha-energy difference is

$$
\Delta Q_{\rm Ds}
=
(11.165-10.858)\frac{273}{269}
=
311.6\ {\rm keV}.
$$

Define

$$
\Delta E_{\rm Ds}
=
E_{\rm Ds}(3/2^+[611])
-
E_{\rm Ds}(13/2^-[716]),
$$

and

$$
\Delta E_{\rm Hs,feed}
=
E_{\rm Hs}
\left[
3/2^+\left(1/2^+[620]\right)
\right]
-
E_{\rm Hs}(11/2^-[725]).
$$

Then

$$
\boxed{
\Delta E_{\rm Ds}
-
\Delta E_{\rm Hs,feed}
=
311.6\ {\rm keV}
}
$$

or equivalently

$$
\boxed{
\Delta E_{\rm Ds}
=
311.6\ {\rm keV}
+
\Delta E_{\rm Hs,feed}.
}
$$

### Drawing consequence

The Ds $3/2^+[611]$ level should be drawn above the $13/2^-[716]$ level.

If the two Hs feeding levels are separated only by several tens of keV, the Ds splitting is naturally of order

$$
\sim 0.3\ {\rm MeV}.
$$

Do not assign an exact Ds separation unless a specific Hs ordering is chosen.

---

# 5. Internal $^{269}$Hs topology

The Hs ordering is allowed to change.

Required low-spin cascade:

$$
3/2^+\left(1/2^+[620]\right)
\rightarrow
1/2^+[620]
\rightarrow
3/2^+[622].
$$

Required high-spin de-excitation:

$$
11/2^-[725]
\rightarrow
9/2^+[615].
$$

Thus only the following local inequalities are mandatory:

$$
E\left[
3/2^+\left(1/2^+[620]\right)
\right]
>
E(1/2^+[620])
>
E(3/2^+[622]),
$$

and

$$
E(11/2^-[725])
>
E(9/2^+[615]).
$$

There is no independent requirement that $9/2^+[615]$ lie above $3/2^+[622]$.

---

# 6. $^{269}$Hs $\rightarrow {}^{265}$Sg: $E_\alpha$ fixes the Hs–Sg relative geometry

Use

$$
^{269}{\rm Hs}(3/2^+[622])
\xrightarrow[\ 9.13\ {\rm MeV}\ ]{\alpha}
^{265}{\rm Sg}_b(3/2^+[622]),
$$

and

$$
^{269}{\rm Hs}(9/2^+[615])
\xrightarrow[\ 8.95\ {\rm MeV}\ ]{\alpha}
^{265}{\rm Sg}_b^{\rm rot}
\left[
9/2^+\left(3/2^+[622]\right)
\right].
$$

The recoil-corrected alpha-energy difference is

$$
\Delta Q_{\rm Hs}
=
(9.13-8.95)\frac{269}{265}
=
182.7\ {\rm keV}.
$$

Define

$$
E_{\rm rot}^{\rm Sg}
=
E_{\rm Sg}
\left[
9/2^+\left(3/2^+[622]\right)
\right]
-
E_{\rm Sg}(3/2^+[622]),
$$

and

$$
\delta_{\rm Hs}
=
E_{\rm Hs}(9/2^+[615])
-
E_{\rm Hs}(3/2^+[622]).
$$

Then

$$
\boxed{
E_{\rm rot}^{\rm Sg}
-
\delta_{\rm Hs}
=
182.7\ {\rm keV}
}
$$

or

$$
\boxed{
E_{\rm rot}^{\rm Sg}
=
182.7\ {\rm keV}
+
\delta_{\rm Hs}.
}
$$

### Drawing consequence

The Hs $9/2^+[615]$ and $3/2^+[622]$ levels cannot be placed arbitrarily once a height is chosen for the Sg rotational member.

Examples only:

- if $E_{\rm rot}^{\rm Sg}=180$ keV, then $\delta_{\rm Hs}\approx-3$ keV;
- if $E_{\rm rot}^{\rm Sg}=150$ keV, then $\delta_{\rm Hs}\approx-33$ keV;
- if $E_{\rm rot}^{\rm Sg}=120$ keV, then $\delta_{\rm Hs}\approx-63$ keV.

Therefore a perfectly acceptable energy-consistent drawing has

$$
E_{\rm Hs}(9/2^+[615])
<
E_{\rm Hs}(3/2^+[622]).
$$

This is why the Hs internal ordering should remain flexible.

---

# 7. $^{265}$Sg $\rightarrow {}^{261}$Rf: long-lived level splittings from $E_\alpha$

Use the long-lived assignments

$$
^{265}{\rm Sg}_a=9/2^+[615],
\qquad
^{265}{\rm Sg}_b=3/2^+[622],
$$

and

$$
^{261}{\rm Rf}_a=9/2^+[615],
\qquad
^{261}{\rm Rf}_b=3/2^+[622].
$$

The observed parent alpha energies are

$$
E_\alpha({\rm Sg}_a)=8.84\ {\rm MeV},
$$

$$
E_\alpha({\rm Sg}_b)=8.69\ {\rm MeV}.
$$

The recoil-corrected difference is

$$
\Delta Q_{\rm Sg}
=
(8.84-8.69)\frac{265}{261}
=
152.3\ {\rm keV}.
$$

For the **diagonal reference assignments**

$$
{\rm Sg}_a\rightarrow{\rm Rf}_a,
$$

and

$$
{\rm Sg}_b\rightarrow{\rm Rf}_b,
$$

define

$$
\Delta E_{\rm Sg}
=
E({\rm Sg}_a)-E({\rm Sg}_b),
$$

and, because Rf$_b$ must be above Rf$_a$,

$$
\Delta E_{\rm Rf}
=
E({\rm Rf}_b)-E({\rm Rf}_a)>0.
$$

Then

$$
\boxed{
\Delta E_{\rm Sg}
+
\Delta E_{\rm Rf}
=
152.3\ {\rm keV}.
}
$$

This is an important energy constraint.

Therefore, if both required orderings hold,

$$
E({\rm Sg}_a)>E({\rm Sg}_b),
$$

and

$$
E({\rm Rf}_b)>E({\rm Rf}_a),
$$

then automatically

$$
0<\Delta E_{\rm Sg}<152.3\ {\rm keV},
$$

and

$$
0<\Delta E_{\rm Rf}<152.3\ {\rm keV}.
$$

### Drawing consequence

Do not draw the Sg$_a$–Sg$_b$ or Rf$_b$–Rf$_a$ separations as several hundred keV.

Their two positive splittings must share the available

$$
152.3\ {\rm keV}.
$$

---

# 8. $^{261}$Rf $\rightarrow {}^{257}$No gives an additional strong constraint

The known No ground state is

$$
^{257}{\rm No}_{\rm g.s.}=3/2^+[622].
$$

Use

$$
^{261}{\rm Rf}_b(3/2^+[622])
\xrightarrow[\ 8.51\ {\rm MeV}\ ]{\alpha}
^{257}{\rm No}_{\rm g.s.}(3/2^+[622]),
$$

and tentatively

$$
^{261}{\rm Rf}_a(9/2^+[615])
\xrightarrow[\ 8.28\ {\rm MeV}\ ]{\alpha}
^{257}{\rm No}^*(9/2^+[615]).
$$

Let

$$
E_{\rm No}^*
=
E_x\left[^{257}{\rm No}(9/2^+[615])\right].
$$

The recoil-corrected alpha-energy difference is

$$
\Delta Q_{\rm Rf}
=
(8.51-8.28)\frac{261}{257}
=
233.6\ {\rm keV}.
$$

Energy conservation gives

$$
\boxed{
\Delta E_{\rm Rf}
+
E_{\rm No}^*
=
233.6\ {\rm keV},
}
$$

so

$$
\boxed{
\Delta E_{\rm Rf}
=
233.6\ {\rm keV}
-
E_{\rm No}^*.
}
$$

Because the required Rf ordering is

$$
E({\rm Rf}_b)>E({\rm Rf}_a),
$$

the tentative No $9/2^+[615]$ state must satisfy

$$
E_{\rm No}^*<233.6\ {\rm keV}.
$$

---

# 9. Combining Sg and Rf constraints: the strongest result for the drawing

From

$$
\Delta E_{\rm Sg}
+
\Delta E_{\rm Rf}
=
152.3\ {\rm keV},
$$

and

$$
\Delta E_{\rm Rf}
+
E_{\rm No}^*
=
233.6\ {\rm keV},
$$

we obtain

$$
\boxed{
\Delta E_{\rm Sg}
=
E_{\rm No}^*
-
81.3\ {\rm keV}.
}
$$

Therefore simultaneous requirements

$$
\Delta E_{\rm Sg}>0
$$

and

$$
\Delta E_{\rm Rf}>0
$$

imply

$$
\boxed{
81.3\ {\rm keV}
<
E_{\rm No}^*
<
233.6\ {\rm keV}.
}
$$

Equivalently,

$$
\boxed{
0<
E({\rm Sg}_a)-E({\rm Sg}_b)
<
152.3\ {\rm keV},
}
$$

$$
\boxed{
0<
E({\rm Rf}_b)-E({\rm Rf}_a)
<
152.3\ {\rm keV},
}
$$

and

$$
\boxed{
\left[E({\rm Sg}_a)-E({\rm Sg}_b)\right]
+
\left[E({\rm Rf}_b)-E({\rm Rf}_a)\right]
=
152.3\ {\rm keV}.
}
$$

This relation should control the relative vertical distances in the Sg and Rf columns.

---

# 10. Recommended schematic coordinates if Codex needs explicit numbers

The data do not determine a unique value of $E_{\rm No}^*$.

If Codex requires explicit coordinates only to draw a readable schematic, a convenient **illustrative, non-physical-fit** choice is

$$
E_{\rm No}^*=150\ {\rm keV}.
$$

Then

$$
E({\rm Rf}_b)-E({\rm Rf}_a)
=
233.6-150
=
83.6\ {\rm keV},
$$

and

$$
E({\rm Sg}_a)-E({\rm Sg}_b)
=
150-81.3
=
68.7\ {\rm keV}.
$$

Thus one energy-consistent illustrative layout is:

### $^{265}$Sg long-lived states

set

$$
E({\rm Sg}_b)=0,
$$

then

$$
E({\rm Sg}_a)\approx69\ {\rm keV}.
$$

### $^{261}$Rf long-lived states

set

$$
E({\rm Rf}_a)=0,
$$

then

$$
E({\rm Rf}_b)\approx84\ {\rm keV}.
$$

### $^{257}$No

$$
E(3/2^+[622])=0,
$$

$$
E(9/2^+[615])\approx150\ {\rm keV}
$$

for this illustrative coordinate choice only.

Do not present $150$, $84$, or $69$ keV as measured energies.

---

# 11. Rotational $9/2^+$ members

The following levels must be shown explicitly:

$$
^{265}{\rm Sg}:
\quad
9/2^+\left(3/2^+[622]\right),
$$

and

$$
^{261}{\rm Rf}:
\quad
9/2^+\left(3/2^+[622]\right).
$$

They must lie above their respective $3/2^+[622]$ band heads:

$$
E\left[
9/2^+\left(3/2^+[622]\right)
\right]
>
E(3/2^+[622]).
$$

For Sg, its height is tied to the Hs splitting through

$$
E_{\rm rot}^{\rm Sg}
=
182.7\ {\rm keV}
+
\delta_{\rm Hs}.
$$

Therefore the previous unconditional ordering

$$
{\rm Sg}_a
>
9/2^+\left(3/2^+[622]\right)
>
{\rm Sg}_b
$$

must **not** be imposed.

Depending on the Hs splitting, the Sg rotational $9/2^+$ member may lie above Sg$_a$.

For Rf, only require

$$
E\left[
9/2^+\left(3/2^+[622]\right)
\right]
>
E({\rm Rf}_b)
>
E({\rm Rf}_a).
$$

---

# 12. The four Sg $\rightarrow$ Rf paths and how Codex should interpret them

The drawing must show all four qualitative relationships:

## Diagonal candidates

$$
{\rm Sg}_a(9/2^+[615])
\rightarrow
{\rm Rf}_a(9/2^+[615]),
$$

$$
{\rm Sg}_b(3/2^+[622])
\rightarrow
{\rm Rf}_b(3/2^+[622]).
$$

## Cross candidate using the rotational member

$$
{\rm Sg}_a(9/2^+[615])
\rightarrow
{\rm Rf}_b^{\rm rot}
\left[
9/2^+\left(3/2^+[622]\right)
\right]
\rightarrow_{\gamma/{\rm IC}}
{\rm Rf}_b(3/2^+[622]).
$$

## Other cross candidate

$$
{\rm Sg}_b(3/2^+[622])
\rightarrow
{\rm Rf}_a(9/2^+[615]).
$$

### Critical energy bookkeeping note

The table uses the same measured Sg$_a$ alpha energy, $8.84$ MeV, when evaluating different possible Rf endpoint assignments, and similarly uses $8.69$ MeV for the Sg$_b$ alternatives.

Therefore the four arrows are **candidate endpoint assignments / competing decay interpretations**.

They must not be interpreted as four independently resolved alpha branches all having those exact energies simultaneously.

If Codex draws all four arrows in one figure, use a legend such as:

- solid / emphasized: preferred diagonal reference assignment;
- dashed: alternative cross assignment;
- $\gamma$/IC cascade after feeding a rotational member.

The long-lived level heights should be determined from the diagonal reference pair, while the cross arrows are overlaid as alternative structural possibilities.

---

# 13. Final qualitative vertical ordering for Codex

## $^{273}$Ds

$$
3/2^+[611]
>
13/2^-[716],
$$

with a splitting of order $\sim0.3$ MeV if the two relevant Hs feeding levels are close.

## $^{269}$Hs

Do not impose one global order.

Require only

$$
3/2^+\left(1/2^+[620]\right)
>
1/2^+[620]
>
3/2^+[622],
$$

and

$$
11/2^-[725]
>
9/2^+[615].
$$

The relative position of $9/2^+[615]$ and $3/2^+[622]$ must be chosen consistently with

$$
E_{\rm rot}^{\rm Sg}
-
\left[
E_{\rm Hs}(9/2^+[615])
-
E_{\rm Hs}(3/2^+[622])
\right]
=
182.7\ {\rm keV}.
$$

## $^{265}$Sg

Mandatory:

$$
E({\rm Sg}_a)>E({\rm Sg}_b),
$$

and

$$
0<E({\rm Sg}_a)-E({\rm Sg}_b)<152.3\ {\rm keV}.
$$

The $9/2^+(3/2^+[622])$ rotational member must be above Sg$_b$, but its position relative to Sg$_a$ is not fixed solely by the present $E_\alpha$ data.

## $^{261}$Rf

Mandatory ordering:

$$
9/2^+\left(3/2^+[622]\right)
>
3/2^+[622]\ ({\rm Rf}_b)
>
9/2^+[615]\ ({\rm Rf}_a).
$$

Also

$$
0<E({\rm Rf}_b)-E({\rm Rf}_a)<152.3\ {\rm keV}.
$$

## $^{257}$No

$$
9/2^+[615]\ {\rm tentative}
>
3/2^+[622]\ {\rm g.s.},
$$

with the energy-consistency window

$$
81.3\ {\rm keV}
<
E_x(9/2^+[615])
<
233.6\ {\rm keV}.
$$

---

# 14. Non-negotiable drawing instructions

1. Do not draw the rejected Cn $\rightarrow$ Ds $l_\alpha=5$ path.
2. Show the two Ds $\rightarrow$ Hs branches.
3. Keep the Hs global ordering flexible.
4. Use the $9.13$ and $8.95$ MeV Hs lines to constrain the Hs/Sg relative heights through the $182.7$ keV relation.
5. Show Sg$_a$, Sg$_b$, and the $9/2^+$ rotational member of the Sg$_b$ band.
6. Show all four Sg $\rightarrow$ Rf candidate relations, including both cross paths.
7. Explicitly show the Rf$_b$-band $9/2^+$ rotational member used by the Sg$_a\rightarrow$Rf$_b$ cross feeding.
8. Draw Rf$_b=3/2^+[622]$ above Rf$_a=9/2^+[615]$.
9. Enforce
   $$
   \Delta E_{\rm Sg}+\Delta E_{\rm Rf}=152.3\ {\rm keV}.
   $$
10. Keep the known $^{257}$No ground state as $3/2^+[622]$.
11. If a tentative No $9/2^+[615]$ level is used to receive Rf$_a$, place it within the energy-consistency interval
   $$
   81.3<E_x<233.6\ {\rm keV}.
   $$
12. Label the whole figure `schematic; relative heights constrained by E_alpha, not an adopted level scheme`.

