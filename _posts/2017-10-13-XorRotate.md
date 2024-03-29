---
layout:       post
title:        Basic XOR-rotates and their inverse
categories:   [math]
tags:         [integer]
description:  Simplified description of bijections formed only from XORs and bit rotations.
---

\\
A quick note on a class of invertible (bijective) integer functions composed solely from XOR and bit-rotation operations. The handful of references I've run across seem to make things way more complicated than needed from my perspective.

* Just assuming 32-bit integers but extends to other widths (but power-of-two only).
* Using convention that a 32-bit integer `x` in matrix form is $x=\left(x_0,x_1,\ldots\right)^{T}$ where $x_0$ is the low bit.

<br>

------

Forward function
------

\\
If you're not familar with the function(s) of this note then the following might be the first thing to come to mind (with `a` some constant):  `x = x ^ rot(x, a)`. 

Sadly this isn't invertible, we can see that by setting `x` to zero and then to -1 and in both cases the result is zero.  A bijection must be a one-to-one mapping.  By adding a second rotate we do get a bijection (`b` is a second constant): `x = x ^ rot(x,a) ^ rot(x,b)`.  Some generic code:

<br>

{% highlight c %}

static inline uint32_t rot(const uint32_t x, uint32_t i)
{
  return (x << i)|(x >> (32-i));
}

static inline uint32_t xor_rot2(uint32_t x, uint32_t a, uint32_t b)
{
  return x^rot(x,a)^rot(x,b);
}

{% endhighlight %}


\\
We can represent this function as a matrix equation over $\mathbb{F}_2$:

$$ x' = M~x$$

\\
where the transform matrix is:

$$ \begin{equation} 
M = I + C^a + C^b  \label{M} 
\end{equation} $$

\\
with $I$ the identity matrix and $C$ the cyclic permutation[^circulant] matrix. This generalizes to a sum of $C$ terms:


$$ \begin{equation} 
M = I + \sum\limits_{i=1}^n C^{k_i} \label{genM}
\end{equation} $$

\\
which is invertible if $n$ is even. Transforms like:

$$ \begin{equation} 
C^a + C^b + C^c \label{rot3}
\end{equation} $$

\\
are bijections as well for odd number of terms, but aren't worth considering (for the moment on the math side) since:


$$ \begin{equation} 
C^a\left(I + C^{b-a} + C^{c-a}\right) \label{rot3f}
\end{equation} $$

\\
We simply have a composition[^composition] of $\eqref{M}$ followed by a bit-rotation.  As code

* $\eqref{rot3}$ `rot(x,a)^rot(x,b)^rot(x,c)` 
* $\eqref{rot3f}$ `rot(x^rot(x,b-a)^rot(x,c-a), a)`

\\
are the same (provided `a` is smaller than `b` and `c`, otherwise simply rename). 

\\
A little CPU wasting eye candy of invertible matrices where the list of numbers indicate the powers of $C$ to be summed:

<p align="center"><canvas id="fig" width="161" height="161" class="center"></canvas>
<div align="center" id="dfig">Foo</div>
</p>

\\
For 32-bit we have 4991 matrices (odd since some are involutions) by generalizing \eqref{genM} and \eqref{rot3} ($i$ still must be even and $k$'s unique):

$$ \begin{equation} 
M = C^{k_0}\left(I + \sum\limits_{i=1}^n C^{k_i}\right)
\end{equation} $$

\\
Some notes on equation manipulation:

* Although working with matrix equations we are limited to a commutative sub-algebra: $M_0~M_1 = M_1~M_0$.
* Powers of $C$ modulo reduce:  $C^{k} = C^{k\bmod 32}$
* Leading coefficients drop since we're working in $\mathbb{F}_2$:  $n~C^k$ = $\left(n \bmod 2\right) C^k$.  Or more simply even become zero and odd become one.

<br>

------

Inverse function
------

\\
Instead of a frontal attack on the inverse function, let's sneak up on it from behind.  First we need to find the period of $M$ which is defined as the smallest positive integer $n$ such that:

$$
M^n = M
$$

\\
Let's start by expanding the function squared in painful detail:

$$ \begin{eqnarray}
M^2 & = & \left(I + C^a + C^b\right)\left(I + C^a + C^b\right)  \label{m2_1} \\
    & = & I^2 + I~C^a + I~C^b + C^a~I + C^a~C^a + C^a~C^b + C^b~I + C^b~C^a +C^b~C^b \label{m2_2} \\ 
    & = & I + C^a + C^b + C^a + C^{2a} + C^{a+b} + C^b + C^{a+b} +C^{2b} \label{m2_3} \\
    & = & I + 2~C^a + 2~C^b +  2~C^{a+b} + C^{2a} + C^{2b} \label{m2_4} \\
    & = & I + C^{2a} + C^{2b} \label{m2} 
\end{eqnarray} $$


\\
with line-by-line commentary:

* just state what's being expanded $\eqref{m2_1}$
* mechanical expansion $\eqref{m2_2}$
* simplify $I$ terms and collect powers of $C$ $\eqref{m2_3}$
* gather common terms $\eqref{m2_4}$
* in $\mathbb{F}_2$ odd coefficients reduce to one and even to zero $\eqref{m2}$

\\
This holds for the generalized version \eqref{genM}. By simply renaming the variables we can repeat this squaring process as much as we like.  So we have:

$$  \begin{eqnarray*}
M^{32} &=& I + C^{32a} + C^{32b} \\
       &=& 3I \\
       &=& I
\end{eqnarray*} $$

\\
By extension we now have a justifiction for the set of expressions which are invertible. Given $n$ terms we have:

$$  \begin{eqnarray*}
M^{32} &=& nI \\
\end{eqnarray*} $$

\\
which reduces to $I$ if $n$ is odd and otherwise to zero and thus $n$ must be odd to be an invertible function. Returning to our example this implies:

$$ \begin{eqnarray}
M^{-1} & = & M^{31} \\
       & = & M~M^{2}~M^{4}~M^{8}~M^{16} \\
       & = & \left(I + C^a + C^b\right)\left(I + C^{2a} + C^{2b}\right)\left(I + C^{4a} + C^{4b}\right)\left(I + C^{8a} + C^{8b}\right)\left(I + C^{16a} + C^{16b}\right) \label{generic}
\end{eqnarray} $$



\\
A generic inverse function (of `xor_rot2` above) is then:

{% highlight c %}
uint32_t xor_rot2_inv(uint32_t x, uint32_t a, uint32_t b)
{
  // Perform the five steps (and keep 'a' and 'b' in range given
  // how 'rot' is defined above) as a sequence of transforms.
  // The order reversed of above (products of powers of M commute).
  x = x^rot(x,a)^rot(x,b); a = (a+a) & 0x1f; b = (b+b) & 0x1f; // t0 = M x
  x = x^rot(x,a)^rot(x,b); a = (a+a) & 0x1f; b = (b+b) & 0x1f; // t1 = M^2 t0
  x = x^rot(x,a)^rot(x,b); a = (a+a) & 0x1f; b = (b+b) & 0x1f; // t2 = M^4 t1
  x = x^rot(x,a)^rot(x,b); a = (a+a) & 0x1f; b = (b+b) & 0x1f; // t3 = M^8 t2
  x = x^rot(x,a)^rot(x,b);                                     // x' = M^16 t3

  return x;
}
{% endhighlight %}

\\
For known constants simply plugging them into $\eqref{generic}$ and modulo reducing the values works but is always performing extra work.  Let's run through an example of finding the inverse of `x^rot(x,11)^rot(x,16)` which gives us:


$$ \begin{eqnarray}
M & = & I + C^{11} + C^{16} \nonumber \\
\nonumber \\
M^{-1} 
& = &  \left(I + C^{11} + C^{16}\right)\left(I + C^{22} + C^{32}\right)\left(I + C^{44} + C^{64}\right)\left(I + C^{88} + C^{128}\right)\left(I + C^{176} + C^{256}\right) \label{x1_1} \\
& = & \left(I + C^{11} + C^{16}\right)\left(I + C^{22} + C^0\right)\left(I + C^{12} + C^0\right)\left(I + C^{24} + C^0\right)\left(I + C^{16} + C^0\right) \label{x1_2} \\
& = & \left(I + C^{11} + C^{16}\right)\left(C^{22}\right)\left(C^{12}\right)\left(C^{24}\right)\left(C^{16}\right) \label{x1_3} \\ 
& = & \left(I + C^{11} + C^{16}\right)\left(C^{10}\right) \label{x1_4} \\
& = & \left(C^{10} + C^{21} + C^{26}\right) \label{x1_5}
\end{eqnarray} $$

\\
Equation $\eqref{x1_1}$ is simply expanding $\eqref{generic}$ and modulo reduction in $\eqref{x1_2}$. Here we have $C^0=I$ which gives us $2I$ in each of the last four product terms.  Recall any even coefficient reduces to zero and any odd to one bring us to $\eqref{x1_3}$. Equation $\eqref{x1_4}$ simply sums up and modulo reduces these terms.  Finally $\eqref{x1_5}$ is the final result which translates back into code as `rot(x,10)^rot(x,21)^rot(x,26)`.

\\
This inverse is in the form of $\eqref{rot3}$ which I dismissed as being uninteresting from a math perspetive.  If this was the function we wanted to inverse then we could factor as per $\eqref{rot3f}$ which would give us $\eqref{x1_4}$, which we'd find the inverse of the left hand side (which ends up being identical to the original) then we'd have to compose that with the inverse of $C^{10}$ which is $C^{-10}$ and back to the original forward transform.



------

Inverse Tables
------

\\
Since I tossed together a Mathematica script to mechanically perform the reductions mentioned above I might as well paste the 465 two rotation equations. This script does not attempt to generate a minimal inverse.  Check any results since I could have screwed-up in the cut/paste/mark-up phase (call it about 5 seconds to exhastively test one pair). 

The tables entries are a list of the powers of $C$ to be summed.  Example: $\left(0,1,17\right) = I+C^1+C^{17}$.

\\
There are 15 self-inverses (involutions, so period is 2):

{: .center }
 | (0, 1, 17)  | (0, 2, 18)  | (0, 3, 19)  | (0, 4, 20) |
 | (0, 5, 21)  | (0, 6, 22)  | (0, 7, 23)  | (0, 8, 24) |
 | (0, 9, 25)  | (0, 10, 26) | (0, 11, 27) |(0, 12, 28) |
 | (0, 13, 29) | (0, 14, 30) | (0, 15, 31) |


\\
Those that reduce to a single product term.  The $n$ is the period of the function.

{: .center }
 |  $M$ | n | $M^{-1}$ |  |  $M$ | n | $M^{-1}$ |
 |----------|-----------|-----------|-----------|-----------|-----------|-----------|
 | (0, 1, 16)  |32 | (14, 30, 31) | | (0, 2, 16)  |16 |(12, 28, 30) |
 | (0, 3, 16)  |32 |(10, 26, 29)  | | (0, 4, 16)  | 8 |(8, 24, 28) |
 | (0, 5, 16)  |32 |(6, 22, 27)   | | (0, 6, 16)  |16 |(4, 20, 26) |
 | (0, 7, 16)  |32 |(2, 18, 25)   | | (0, 8, 16)  | 4 |(0, 16, 24) |
 | (0, 9, 16)  |32 |(14, 23, 30)  | | (0, 10, 16) |16 |(12, 22, 28) |
 | (0, 11, 16) |32 |(10, 21, 26)  | | (0, 12, 16) | 8 |(8, 20, 24) |
 | (0, 13, 16) |32 |(6, 19, 22)   | | (0, 14, 16) |16 |(4, 18, 20) |
 | (0, 15, 16) |32 |(2, 17, 18)   | | (0, 16, 17) |32 |(14, 15, 30) |
 | (0, 16, 18) |16 |(12, 14, 28)  | | (0, 16, 19) |32 |(10, 13, 26) |
 | (0, 16, 20) | 8 |(8, 12, 24)   | | (0, 16, 21) |32 |(6, 10, 22) |
 | (0, 16, 22) |16 |(4, 10, 20)   | | (0, 16, 23) |32 |(2, 9, 18) |
 | (0, 16, 24) | 4 |(0, 8, 16)    | | (0, 16, 25) |32 |(7, 14, 30) |
 | (0, 16, 26) |16 |(6, 12, 28)   | | (0, 16, 27) |32 |(5, 10, 26) |
 | (0, 16, 28) | 8 |(4, 8, 24)    | | (0, 16, 29) |32 |(3, 6, 22) |
 | (0, 16, 30) |16 |(2, 4, 20)    | | (0, 16, 31) |32 |(1, 2, 18) |

\\
The two product term inverses:

{: .center }
 |  $M$ | n | $M^{-1}$ |  |  $M$ | n | $M^{-1}$ |
 |----------|-----------|-----------|-----------|-----------|-----------|-----------|
 | (0, 1, 8) | 32| (0, 1, 8)(12, 28, 30) | | (0, 2, 8) | 16| (0, 4, 16)(0, 24, 26) |
 | (0, 3, 8) | 32| (0, 3, 8)(4, 20, 26) | | (0, 4, 8) | 8| (0, 4, 8)(0, 16, 24) |
 | (0, 5, 8) | 32| (0, 5, 8)(12, 22, 28) || (0, 6, 8) | 16| (0, 6, 8)(8, 20, 24) |
 | (0, 7, 8) | 32| (0, 7, 8)(4, 18, 20) | | (0, 1, 9) | 4| (0, 1, 9)(0, 2, 18) |
 | (0, 8, 9) | 32| (0, 8, 9)(12, 14, 28) | | (0, 2, 10) | 4| (0, 2, 10)(0, 4, 20) |
 | (0, 8, 10) | 16| (0, 16, 20)(0, 2, 24) | | (0, 3, 11) | 4| (0, 3, 11)(0, 6, 22) |
 | (0, 8, 11) | 32| (0, 8, 11)(4, 10, 20) | | (0, 4, 12) | 4| (0, 4, 12)(0, 8, 24) |
 | (0, 8, 12) | 8| (0, 8, 12)(0, 8, 16) | | (0, 5, 13) | 4| (0, 5, 13)(0, 10, 26) |
 | (0, 8, 13) | 32| (0, 8, 13)(6, 12, 28) | | (0, 6, 14) | 4| (0, 6, 14)(0, 12, 28) |
 | (0, 8, 14) | 16| (0, 8, 14)(4, 8, 24) | | (0, 7, 15) | 4| (0, 7, 15)(0, 14, 30) |
 | (0, 8, 15) | 32| (0, 8, 15)(2, 4, 20) | | (0, 8, 17) | 32| (0, 2, 16)(4, 13, 28) |
 | (0, 9, 17) | 4| (0, 9, 17)(0, 2, 18) | | (0, 8, 18) | 16| (0, 4, 16)(0, 10, 24) |
 | (0, 10, 18) | 4| (0, 10, 18)(0, 4, 20) | | (0, 8, 19) | 32| (0, 6, 16)(7, 20, 28) |
 | (0, 11, 19) | 4| (0, 11, 19)(0, 6, 22) | | (0, 8, 20) | 8| (0, 8, 20)(0, 16, 24) |
 | (0, 12, 20) | 4| (0, 12, 20)(0, 8, 24) | | (0, 8, 21) | 32| (0, 10, 16)(1, 12, 20) |
 | (0, 13, 21) | 4| (0, 13, 21)(0, 10, 26) | | (0, 8, 22) | 16| (0, 12, 16)(8, 16, 30) |
 | (0, 14, 22) | 4| (0, 14, 22)(0, 12, 28) | | (0, 8, 23) | 32| (0, 14, 16)(4, 12, 27) |
 | (0, 15, 23) | 4| (0, 15, 23)(0, 14, 30) | | (0, 1, 24) | 32| (0, 2, 16)(20, 28, 29) |
 | (0, 2, 24) | 16| (0, 4, 16)(16, 24, 26) | | (0, 3, 24) | 32| (0, 6, 16)(12, 20, 23) |
 | (0, 4, 24) | 8| (0, 4, 24)(0, 16, 24) | | (0, 5, 24) | 32| (0, 10, 16)(4, 12, 17) |
 | (0, 6, 24) | 16| (0, 8, 14)(0, 12, 16) | | (0, 7, 24) | 32| (0, 14, 16)(4, 11, 28) |
 | (0, 9, 24) | 32| (0, 16, 18)(5, 20, 28) | | (0, 10, 24) | 16| (0, 16, 20)(2, 16, 24) |
 | (0, 11, 24) | 32| (0, 16, 22)(12, 20, 31) | | (0, 12, 24) | 8| (0, 8, 16)(0, 12, 24) |
 | (0, 13, 24) | 32| (0, 13, 24)(6, 12, 28) | | (0, 14, 24) | 16| (0, 8, 22)(0, 16, 28) |
 | (0, 15, 24) | 32| (2, 4, 20)(0, 15, 24) | | (0, 17, 24) | 32| (0, 2, 16)(13, 20, 28) |
 | (0, 18, 24) | 16| (0, 4, 16)(10, 16, 24) | | (0, 19, 24) | 32| (0, 6, 16)(7, 12, 20) |
 | (0, 20, 24) | 8| (0, 16, 24)(0, 20, 24) | | (0, 21, 24) | 32| (1, 4, 12)(0, 10, 16) |
 | (0, 22, 24) | 16| (0, 12, 16)(0, 8, 30) | | (0, 23, 24) | 32| (0, 14, 16)(4, 27, 28) |
 | (0, 1, 25) | 4| (0, 2, 18)(0, 1, 25) | | (0, 8, 25) | 32| (0, 16, 18)(4, 21, 28) |
 | (0, 17, 25) | 4| (0, 2, 18)(0, 17, 25) | | (0, 24, 25) | 32| (0, 16, 18)(20, 21, 28) |
 | (0, 2, 26) | 4| (0, 4, 20)(0, 2, 26) | | (0, 8, 26) | 16| (0, 16, 20)(0, 18, 24) |
 | (0, 18, 26) | 4| (0, 4, 20)(0, 18, 26) | | (0, 24, 26) | 16| (0, 16, 20)(16, 18, 24) |
 | (0, 3, 27) | 4| (0, 6, 22)(0, 3, 27) | | (0, 8, 27) | 32| (0, 16, 22)(15, 20, 28) |
 | (0, 19, 27) | 4| (0, 6, 22)(0, 19, 27) | | (0, 24, 27) | 32| (12, 15, 20)(0, 16, 22) |
 | (0, 4, 28) | 4| (0, 8, 24)(0, 4, 28) | | (0, 8, 28) | 8| (0, 8, 16)(0, 8, 28) |
 | (0, 20, 28) | 4| (0, 8, 24)(0, 20, 28) | | (0, 24, 28) | 8| (0, 8, 16)(0, 24, 28) |
 | (0, 5, 29) | 4| (0, 10, 26)(0, 5, 29) | | (0, 8, 29) | 32| (9, 12, 20)(0, 16, 26) |
 | (0, 21, 29) | 4| (0, 10, 26)(0, 21, 29) | | (0, 24, 29) | 32| (4, 9, 12)(0, 16, 26) |
 | (0, 6, 30) | 4| (0, 12, 28)(0, 6, 30) | | (0, 8, 30) | 16| (6, 8, 16)(0, 16, 28) |
 | (0, 22, 30) | 4| (0, 12, 28)(0, 22, 30) | | (0, 24, 30) | 16| (0, 6, 8)(0, 16, 28) |
 | (0, 7, 31) | 4| (0, 14, 30)(0, 7, 31) | | (0, 8, 31) | 32| (3, 4, 12)(0, 16, 30) |
 | (0, 23, 31) | 4| (0, 14, 30)(0, 23, 31) | | (0, 24, 31) | 32| (3, 4, 28)(0, 16, 30) |

\\
The three product term inverses:

{: .center }
 |  $M$ | n | $M^{-1}$ |
 |---------|-----------|-----------|
 | (0, 1, 4) | 32| (0, 1, 4)(0, 4, 16)(0, 24, 26) |
 | (0, 2, 4) | 16| (0, 2, 4)(0, 4, 8)(0, 16, 24) |
 | (0, 3, 4) | 32| (0, 3, 4)(0, 6, 8)(8, 20, 24) |
 | (0, 1, 5) | 8| (0, 1, 5)(0, 2, 10)(0, 4, 20) |
 | (0, 4, 5) | 32| (0, 4, 5)(0, 16, 20)(0, 2, 24) |
 | (0, 2, 6) | 8| (0, 2, 6)(0, 4, 12)(0, 8, 24) |
 | (0, 4, 6) | 16| (0, 4, 6)(0, 8, 12)(0, 8, 16) |
 | (0, 3, 7) | 8| (0, 3, 7)(0, 6, 14)(0, 12, 28) |
 | (0, 4, 7) | 32| (0, 4, 7)(0, 8, 14)(4, 8, 24) |
 | (0, 4, 9) | 32| (0, 4, 9)(0, 4, 16)(0, 10, 24) |
 | (0, 5, 9) | 8| (0, 5, 9)(0, 10, 18)(0, 4, 20) |
 | (0, 4, 10) | 16| (0, 4, 10)(0, 8, 20)(0, 16, 24) |
 | (0, 6, 10) | 8| (0, 6, 10)(0, 12, 20)(0, 8, 24) |
 | (0, 4, 11) |  32| (0, 4, 11)(0, 12, 16)(8, 16, 30) |
 | (0, 7, 11) | 8| (0, 7, 11)(0, 14, 22)(0, 12, 28) |
 | (0, 1, 12) |  32| (0, 1, 12)(0, 4, 16)(16, 24, 26) |
 | (0, 2, 12) | 16| (0, 2, 12)(0, 4, 24)(0, 16, 24) |
 | (0, 3, 12) | 32| (0, 3, 12)(0, 8, 14)(0, 12, 16) |
 | (0, 5, 12) |  32| (0, 5, 12)(0, 16, 20)(2, 16, 24) |
 | (0, 6, 12) | 16| (0, 6, 12)(0, 8, 16)(0, 12, 24) |
 | (0, 7, 12) | 32| (0, 7, 12)(0, 8, 22)(0, 16, 28) |
 | (0, 9, 12) |  32| (0, 9, 12)(0, 4, 16)(10, 16, 24) |
 | (0, 10, 12) |  16| (0, 10, 12)(0, 16, 24)(0, 20, 24) |
 | (0, 11, 12) |  32| (0, 11, 12)(0, 12, 16)(0, 8, 30) |
 | (0, 1, 13) | 8| (0, 1, 13)(0, 4, 20)(0, 2, 26) |
 | (0, 4, 13) |  32| (0, 4, 13)(0, 16, 20)(0, 18, 24) |
 | (0, 9, 13) | 8| (0, 9, 13)(0, 4, 20)(0, 18, 26) |
 | (0, 12, 13) |  32| (0, 12, 13)(0, 16, 20)(16, 18, 24) |
 | (0, 2, 14) | 8| (0, 2, 14)(0, 8, 24)(0, 4, 28) |
 | (0, 4, 14) | 16| (0, 4, 14)(0, 8, 16)(0, 8, 28) |
 | (0, 10, 14) |  8| (0, 10, 14)(0, 8, 24)(0, 20, 28) |
 | (0, 12, 14) |  16| (0, 12, 14)(0, 8, 16)(0, 24, 28) |
 | (0, 3, 15) | 8| (0, 3, 15)(0, 12, 28)(0, 6, 30) |
 | (0, 4, 15) |  32| (0, 4, 15)(6, 8, 16)(0, 16, 28) |
 | (0, 11, 15) |  8| (0, 11, 15)(0, 12, 28)(0, 22, 30) |
 | (0, 12, 15) |  32| (0, 6, 8)(0, 12, 15)(0, 16, 28) |
 | (0, 4, 17) | 32| (0, 4, 16)(0, 4, 17)(0, 24, 26) |
 | (0, 5, 17) | 8| (0, 2, 10)(0, 5, 17)(0, 4, 20) |
 | (0, 12, 17) |  32| (0, 4, 16)(0, 12, 17)(16, 24, 26) |
 | (0, 13, 17) | 8| (0, 13, 17)(0, 4, 20)(0, 2, 26) |
 | (0, 4, 18) | 16| (0, 4, 8)(0, 4, 18)(0, 16, 24) |
 | (0, 6, 18) | 8| (0, 4, 12)(0, 6, 18)(0, 8, 24) |
 | (0, 12, 18) |  16| (0, 12, 18)(0, 4, 24)(0, 16, 24) |
 | (0, 14, 18) | 8| (0, 14, 18)(0, 8, 24)(0, 4, 28) |
 | (0, 4, 19) |  32| (0, 6, 8)(0, 12, 16)(8, 12, 27) |
 | (0, 7, 19) | 8| (0, 6, 14)(0, 7, 19)(0, 12, 28) |
 | (0, 12, 19) |  32| (0, 8, 14)(0, 12, 16)(0, 12, 19) |
 | (0, 15, 19) | 8| (0, 15, 19)(0, 12, 28)(0, 6, 30) |
 | (0, 1, 20) | 32| (0, 4, 16)(0, 1, 20)(0, 24, 26) |
 | (0, 2, 20) | 16| (0, 4, 8)(0, 2, 20)(0, 16, 24) |
 | (0, 3, 20) |  32| (0, 6, 8)(0, 12, 16)(8, 11, 28) |
 | (0, 5, 20) | 32| (0, 5, 20)(0, 16, 20)(0, 2, 24) |
 | (0, 6, 20) | 16| (0, 8, 12)(0, 8, 16)(0, 6, 20) |
 | (0, 7, 20) |  32| (0, 8, 14)(0, 7, 20)(4, 8, 24) |
 | (0, 9, 20) | 32| (0, 4, 16)(0, 9, 20)(0, 10, 24) |
 | (0, 10, 20) |  16| (0, 8, 20)(0, 10, 20)(0, 16, 24) |
 | (0, 11, 20) |  32| (0, 12, 16)(0, 11, 20)(8, 16, 30) |
 | (0, 13, 20) |  32| (0, 13, 20)(0, 16, 20)(0, 18, 24) |
 | (0, 14, 20) |  16| (0, 8, 16)(0, 14, 20)(0, 8, 28) |
 | (0, 15, 20) |  32| (6, 8, 16)(0, 15, 20)(0, 16, 28) |
 | (0, 17, 20) |  32| (0, 4, 16)(0, 17, 20)(0, 24, 26) |
 | (0, 18, 20) |  16| (0, 4, 8)(0, 18, 20)(0, 16, 24) |
 | (0, 19, 20) |  32| (0, 6, 8)(0, 12, 16)(8, 27, 28) |
 | (0, 1, 21) | 8| (0, 2, 10)(0, 4, 20)(0, 1, 21) |
 | (0, 4, 21) | 32| (0, 16, 20)(0, 4, 21)(0, 2, 24) |
 | (0, 9, 21) | 8| (0, 10, 18)(0, 4, 20)(0, 9, 21) |
 | (0, 12, 21) |  32| (0, 16, 20)(0, 12, 21)(2, 16, 24) |
 | (0, 17, 21) | 8| (0, 2, 10)(0, 4, 20)(0, 17, 21) |
 | (0, 20, 21) |  32| (0, 16, 20)(0, 20, 21)(0, 2, 24) |
 | (0, 2, 22) | 8| (0, 4, 12)(0, 2, 22)(0, 8, 24) |
 | (0, 4, 22) | 16| (0, 8, 12)(0, 8, 16)(0, 4, 22) |
 | (0, 10, 22) |  8| (0, 12, 20)(0, 10, 22)(0, 8, 24) |
 | (0, 12, 22) |  16| (0, 8, 16)(0, 12, 22)(0, 12, 24) |
 | (0, 18, 22) | 8| (0, 4, 12)(0, 18, 22)(0, 8, 24) |
 | (0, 20, 22) |  16| (0, 8, 12)(0, 8, 16)(0, 20, 22) |
 | (0, 3, 23) | 8| (0, 6, 14)(0, 3, 23)(0, 12, 28) |
 | (0, 4, 23) |  32| (0, 8, 14)(0, 4, 23)(4, 8, 24) |
 | (0, 11, 23) |  8| (0, 14, 22)(0, 11, 23)(0, 12, 28) |
 | (0, 12, 23) |  32| (0, 8, 22)(0, 12, 23)(0, 16, 28) |
 | (0, 19, 23) |  8| (0, 6, 14)(0, 19, 23)(0, 12, 28) |
 | (0, 20, 23) |  32| (0, 8, 14)(0, 20, 23)(4, 8, 24) |
 | (0, 4, 25) | 32| (0, 4, 16)(0, 10, 24)(0, 4, 25) |
 | (0, 5, 25) | 8| (0, 10, 18)(0, 4, 20)(0, 5, 25) |
 | (0, 12, 25) |  32| (0, 4, 16)(4, 17, 24)(0, 18, 24) |
 | (0, 13, 25) |  8| (0, 4, 20)(0, 13, 25)(0, 18, 26) |
 | (0, 20, 25) |  32| (0, 4, 16)(0, 10, 24)(0, 20, 25) |
 | (0, 21, 25) |  8| (0, 10, 18)(0, 4, 20)(0, 21, 25) |
 | (0, 4, 26) | 16| (0, 8, 20)(0, 16, 24)(0, 4, 26) |
 | (0, 6, 26) | 8| (0, 12, 20)(0, 8, 24)(0, 6, 26) |
 | (0, 12, 26) |  16| (0, 16, 24)(0, 20, 24)(0, 12, 26) |
 | (0, 14, 26) |  8| (0, 8, 24)(0, 14, 26)(0, 20, 28) |
 | (0, 20, 26) |  16| (0, 8, 20)(0, 16, 24)(0, 20, 26) |
 | (0, 22, 26) |  8| (0, 12, 20)(0, 8, 24)(0, 22, 26) |
 | (0, 4, 27) |  32| (3, 8, 12)(0, 12, 16)(0, 8, 22) |
 | (0, 7, 27) | 8| (0, 14, 22)(0, 7, 27)(0, 12, 28) |
 | (0, 12, 27) |  32| (0, 12, 16)(0, 12, 27)(0, 8, 30) |
 | (0, 15, 27) |  8| (0, 15, 27)(0, 12, 28)(0, 22, 30) |
 | (0, 20, 27) | 32| (0, 12, 16)(0, 8, 22)(3, 8, 28) |
 | (0, 23, 27) |  8| (0, 14, 22)(0, 23, 27)(0, 12, 28) |
 | (0, 1, 28) |  32| (0, 4, 16)(0, 2, 24)(20, 24, 25) |
 | (0, 2, 28) | 16| (0, 4, 24)(0, 16, 24)(0, 2, 28) |
 | (0, 3, 28) | 32| (0, 8, 14)(0, 12, 16)(0, 3, 28) |
 | (0, 5, 28) |  32| (0, 16, 20)(0, 10, 24)(20, 24, 29) |
 | (0, 6, 28) | 16| (0, 8, 16)(0, 12, 24)(0, 6, 28) |
 | (0, 7, 28) | 32| (0, 8, 22)(0, 7, 28)(0, 16, 28) |
 | (0, 9, 28) |  32| (0, 4, 16)(0, 18, 24)(1, 20, 24) |
 | (0, 10, 28) |  16| (0, 16, 24)(0, 20, 24)(0, 10, 28) |
 | (0, 11, 28) |  32| (0, 12, 16)(0, 11, 28)(0, 8, 30) |
 | (0, 13, 28) |  32| (0, 16, 20)(5, 20, 24)(0, 24, 26) |
 | (0, 14, 28) |  16| (0, 8, 16)(0, 14, 28)(0, 24, 28) |
 | (0, 15, 28) |  32| (0, 6, 8)(0, 15, 28)(0, 16, 28) |
 | (0, 17, 28) |  32| (0, 4, 16)(0, 2, 24)(9, 20, 24) |
 | (0, 18, 28) |  16| (0, 4, 24)(0, 16, 24)(0, 18, 28) |
 | (0, 19, 28) |  32| (0, 8, 14)(0, 12, 16)(0, 19, 28) |
 | (0, 21, 28) |  32| (0, 16, 20)(0, 10, 24)(13, 20, 24) |
 | (0, 22, 28) |  16| (0, 8, 16)(0, 12, 24)(0, 22, 28) |
 | (0, 23, 28) |  32| (0, 8, 22)(0, 16, 28)(0, 23, 28) |
 | (0, 25, 28) |  32| (0, 4, 16)(0, 18, 24)(17, 20, 24) |
 | (0, 26, 28) |  16| (0, 16, 24)(0, 20, 24)(0, 26, 28) |
 | (0, 27, 28) |  32| (0, 12, 16)(0, 27, 28)(0, 8, 30) |
 | (0, 1, 29) | 8| (0, 4, 20)(0, 2, 26)(0, 1, 29) |
 | (0, 4, 29) |  32| (0, 16, 20)(0, 18, 24)(0, 4, 29) |
 | (0, 9, 29) | 8| (0, 4, 20)(0, 18, 26)(0, 9, 29) |
 | (0, 12, 29) |  32| (0, 16, 20)(4, 21, 24)(0, 24, 26) |
 | (0, 17, 29) | 8| (0, 4, 20)(0, 2, 26)(0, 17, 29) |
 | (0, 20, 29) |  32| (0, 16, 20)(0, 18, 24)(0, 20, 29) |
 | (0, 25, 29) |  8| (0, 4, 20)(0, 18, 26)(0, 25, 29) |
 | (0, 28, 29) |  32| (0, 16, 20)(20, 21, 24)(0, 24, 26) |
 | (0, 2, 30) | 8| (0, 8, 24)(0, 4, 28)(0, 2, 30) |
 | (0, 4, 30) | 16| (0, 8, 16)(0, 8, 28)(0, 4, 30) |
 | (0, 10, 30) |  8| (0, 8, 24)(0, 20, 28)(0, 10, 30) |
 | (0, 12, 30) |  16| (0, 8, 16)(0, 24, 28)(0, 12, 30) |
 | (0, 18, 30) | 8| (0, 8, 24)(0, 4, 28)(0, 18, 30) |
 | (0, 20, 30) |  16| (0, 8, 16)(0, 8, 28)(0, 20, 30) |
 | (0, 26, 30) |  8| (0, 8, 24)(0, 20, 28)(0, 26, 30) |
 | (0, 28, 30) |  16| (0, 8, 16)(0, 24, 28)(0, 28, 30) |
 | (0, 3, 31) | 8| (0, 12, 28)(0, 6, 30)(0, 3, 31) |
 | (0, 4, 31) |  32| (7, 8, 12)(0, 16, 28)(0, 8, 30) |
 | (0, 11, 31) |  8| (0, 12, 28)(0, 22, 30)(0, 11, 31) |
 | (0, 12, 31) |  32| (0, 6, 8)(0, 16, 28)(0, 12, 31) |
 | (0, 19, 31) |  8| (0, 12, 28)(0, 6, 30)(0, 19, 31) |
 | (0, 20, 31) |  32| (7, 8, 28)(0, 16, 28)(0, 8, 30) |
 | (0, 27, 31) |  8| (0, 12, 28)(0, 22, 30)(0, 27, 31) |
 | (0, 28, 31) | 32| (0, 6, 8)(0, 16, 28)(0, 28, 31) |


\\
and finally the four product term inverses:

{: .center }
 |  $M$ | n | $M^{-1}$ |
 |------------|-----------|-----------|
 | (0,  1, 2) | 32| (0, 1, 2)(0, 2, 4)(0, 4, 8)(0, 16, 24) |
 | (0,  1, 3) | 16| (0, 1, 3)(0, 2, 6)(0, 4, 12)(0, 8, 24) |
 | (0,  2, 3) | 32| (0, 2, 3)(0, 4, 6)(0, 8, 12)(0, 8, 16) |
 | (0,  2, 5) | 32| (0, 2, 5)(0, 4, 10)(0, 8, 20)(0, 16, 24) |
 | (0,  3, 5) | 16| (0, 3, 5)(0, 6, 10)(0, 12, 20)(0, 8, 24) |
 | (0,  1, 6) | 32| (0, 1, 6)(0, 2, 12)(0, 4, 24)(0, 16, 24) |
 | (0,  3, 6) | 32| (0, 3, 6)(0, 6, 12)(0, 8, 16)(0, 12, 24) |
 | (0,  5, 6) | 32| (0, 5, 6)(0, 10, 12)(0, 16, 24)(0, 20, 24) |
 | (0,  1, 7) | 16| (0, 1, 7)(0, 2, 14)(0, 8, 24)(0, 4, 28) |
 | (0,  2, 7) | 32| (0, 2, 7)(0, 4, 14)(0, 8, 16)(0, 8, 28) |
 | (0,  5, 7) | 16| (0, 5, 7)(0, 10, 14)(0, 8, 24)(0, 20, 28) |
 | (0,  6, 7) | 32| (0, 6, 7)(0, 12, 14)(0, 8, 16)(0, 24, 28) |
 | (0,  2, 9) | 32| (0, 4, 8)(0, 2, 9)(0, 4, 18)(0, 16, 24) |
 | (0,  3, 9) | 16| (0, 3, 9)(0, 4, 12)(0, 6, 18)(0, 8, 24) |
 | (0,  6, 9) | 32| (0, 6, 9)(0, 12, 18)(0, 4, 24)(0, 16, 24) |
 | (0,  7, 9) | 16| (0, 7, 9)(0, 14, 18)(0, 8, 24)(0, 4, 28) |
 | (0,  1, 10) | 32| (0, 4, 8)(0, 1, 10)(0, 2, 20)(0, 16, 24) |
 | (0,  3, 10) | 32| (0, 3, 10)(0, 8, 12)(0, 8, 16)(0, 6, 20) |
 | (0,  5, 10) | 32| (0, 5, 10)(0, 8, 20)(0, 10, 20)(0, 16, 24) |
 | (0,  7, 10) | 32| (0, 7, 10)(0, 8, 16)(0, 14, 20)(0, 8, 28) |
 | (0,  9, 10) | 32| (0, 4, 8)(0, 9, 10)(0, 18, 20)(0, 16, 24) |
 | (0,  1, 11) | 16| (0, 1, 11)(0, 4, 12)(0, 2, 22)(0, 8, 24) |
 | (0,  2, 11) | 32| (0, 2, 11)(0, 8, 12)(0, 8, 16)(0, 4, 22) |
 | (0,  5, 11) | 16| (0, 5, 11)(0, 12, 20)(0, 10, 22)(0, 8, 24) |
 | (0,  6, 11) | 32| (0, 6, 11)(0, 8, 16)(0, 12, 22)(0, 12, 24) |
 | (0,  9, 11) | 16| (0, 9, 11)(0, 4, 12)(0, 18, 22)(0, 8, 24) |
 | (0, 10, 11) | 32| (0, 10, 11)(0, 8, 12)(0, 8, 16)(0, 20, 22) |
 | (0,  2, 13) | 32| (0, 2, 13)(0, 8, 20)(0, 16, 24)(0, 4, 26) |
 | (0,  3, 13) | 16| (0, 3, 13)(0, 12, 20)(0, 8, 24)(0, 6, 26) |
 | (0,  6, 13) | 32| (0, 6, 13)(0, 16, 24)(0, 20, 24)(0, 12, 26) |
 | (0,  7, 13) | 16| (0, 7, 13)(0, 8, 24)(0, 14, 26)(0, 20, 28) |
 | (0, 10, 13) | 32| (0, 10, 13)(0, 8, 20)(0, 16, 24)(0, 20, 26) |
 | (0, 11, 13) | 16| (0, 11, 13)(0, 12, 20)(0, 8, 24)(0, 22, 26) |
 | (0,  1, 14) | 32| (0, 1, 14)(0, 4, 24)(0, 16, 24)(0, 2, 28) |
 | (0,  3, 14) | 32| (0, 3, 14)(0, 8, 16)(0, 12, 24)(0, 6, 28) |
 | (0,  5, 14) | 32| (0, 5, 14)(0, 16, 24)(0, 20, 24)(0, 10, 28) |
 | (0,  7, 14) | 32| (0, 7, 14)(0, 8, 16)(0, 14, 28)(0, 24, 28) |
 | (0,  9, 14) | 32| (0, 9, 14)(0, 4, 24)(0, 16, 24)(0, 18, 28) |
 | (0, 11, 14) | 32| (0, 11, 14)(0, 8, 16)(0, 12, 24)(0, 22, 28) |
 | (0, 13, 14) | 32| (0, 13, 14)(0, 16, 24)(0, 20, 24)(0, 26, 28) |
 | (0,  1, 15) | 16| (0, 1, 15)(0, 8, 24)(0, 4, 28)(0, 2, 30) |
 | (0,  2, 15) | 32| (0, 2, 15)(0, 8, 16)(0, 8, 28)(0, 4, 30) |
 | (0,  5, 15) | 16| (0, 5, 15)(0, 8, 24)(0, 20, 28)(0, 10, 30) |
 | (0,  6, 15) | 32| (0, 6, 15)(0, 8, 16)(0, 24, 28)(0, 12, 30) |
 | (0,  9, 15) | 16| (0, 9, 15)(0, 8, 24)(0, 4, 28)(0, 18, 30) |
 | (0, 10, 15) | 32| (0, 10, 15)(0, 8, 16)(0, 8, 28)(0, 20, 30) |
 | (0, 13, 15) | 16| (0, 13, 15)(0, 8, 24)(0, 20, 28)(0, 26, 30) |
 | (0, 14, 15) | 32| (0, 14, 15)(0, 8, 16)(0, 24, 28)(0, 28, 30) |
 | (0,  2, 17) | 32| (0, 2, 4)(0, 4, 8)(0, 2, 17)(0, 16, 24) |
 | (0,  3, 17) | 16| (0, 2, 6)(0, 4, 12)(0, 3, 17)(0, 8, 24) |
 | (0,  6, 17) | 32| (0, 2, 12)(0, 6, 17)(0, 4, 24)(0, 16, 24) |
 | (0,  7, 17) | 16| (0, 2, 14)(0, 7, 17)(0, 8, 24)(0, 4, 28) |
 | (0, 10, 17) | 32| (0, 4, 8)(0, 10, 17)(0, 2, 20)(0, 16, 24) |
 | (0, 11, 17) | 16| (0, 4, 12)(0, 11, 17)(0, 2, 22)(0, 8, 24) |
 | (0, 14, 17) | 32| (0, 14, 17)(0, 4, 24)(0, 16, 24)(0, 2, 28) |
 | (0, 15, 17) | 16| (0, 15, 17)(0, 8, 24)(0, 4, 28)(0, 2, 30) |
 | (0,  1, 18) | 32| (0, 2, 4)(0, 4, 8)(0, 1, 18)(0, 16, 24) |
 | (0,  3, 18) | 32| (0, 4, 6)(0, 8, 12)(0, 8, 16)(0, 3, 18) |
 | (0,  5, 18) | 32| (0, 4, 10)(0, 5, 18)(0, 8, 20)(0, 16, 24) |
 | (0,  7, 18) | 32| (0, 4, 14)(0, 8, 16)(0, 7, 18)(0, 8, 28) |
 | (0,  9, 18) | 32| (0, 4, 8)(0, 4, 18)(0, 9, 18)(0, 16, 24) |
 | (0, 11, 18) | 32| (0, 8, 12)(0, 8, 16)(0, 11, 18)(0, 4, 22) |
 | (0, 13, 18) | 32| (0, 13, 18)(0, 8, 20)(0, 16, 24)(0, 4, 26) |
 | (0, 15, 18) | 32| (0, 8, 16)(0, 15, 18)(0, 8, 28)(0, 4, 30) |
 | (0, 17, 18) | 32| (0, 2, 4)(0, 4, 8)(0, 17, 18)(0, 16, 24) |
 | (0,  1, 19) | 16| (0, 2, 6)(0, 4, 12)(0, 1, 19)(0, 8, 24) |
 | (0,  2, 19) | 32| (0, 4, 6)(0, 8, 12)(0, 8, 16)(0, 2, 19) |
 | (0,  5, 19) | 16| (0, 6, 10)(0, 5, 19)(0, 12, 20)(0, 8, 24) |
 | (0,  6, 19) | 32| (0, 6, 12)(0, 8, 16)(0, 6, 19)(0, 12, 24) |
 | (0,  9, 19) | 16| (0, 4, 12)(0, 6, 18)(0, 9, 19)(0, 8, 24) |
 | (0, 10, 19) | 32| (0, 8, 12)(0, 8, 16)(0, 10, 19)(0, 6, 20) |
 | (0, 13, 19) | 16| (0, 13, 19)(0, 12, 20)(0, 8, 24)(0, 6, 26) |
 | (0, 14, 19) | 32| (0, 8, 16)(0, 14, 19)(0, 12, 24)(0, 6, 28) |
 | (0, 17, 19) | 16| (0, 2, 6)(0, 4, 12)(0, 17, 19)(0, 8, 24) |
 | (0, 18, 19) | 32| (0, 4, 6)(0, 8, 12)(0, 8, 16)(0, 18, 19) |
 | (0,  2, 21) | 32| (0, 4, 10)(0, 8, 20)(0, 2, 21)(0, 16, 24) |
 | (0,  3, 21) | 16| (0, 6, 10)(0, 12, 20)(0, 3, 21)(0, 8, 24) |
 | (0,  6, 21) | 32| (0, 10, 12)(0, 6, 21)(0, 16, 24)(0, 20, 24) |
 | (0,  7, 21) | 16| (0, 10, 14)(0, 7, 21)(0, 8, 24)(0, 20, 28) |
 | (0, 10, 21) | 32| (0, 8, 20)(0, 10, 20)(0, 10, 21)(0, 16, 24) |
 | (0, 11, 21) | 16| (0, 12, 20)(0, 11, 21)(0, 10, 22)(0, 8, 24) |
 | (0, 14, 21) | 32| (0, 14, 21)(0, 16, 24)(0, 20, 24)(0, 10, 28) |
 | (0, 15, 21) | 16| (0, 15, 21)(0, 8, 24)(0, 20, 28)(0, 10, 30) |
 | (0, 18, 21) | 32| (0, 4, 10)(0, 8, 20)(0, 18, 21)(0, 16, 24) |
 | (0, 19, 21) | 16| (0, 6, 10)(0, 12, 20)(0, 19, 21)(0, 8, 24) |
 | (0,  1, 22) | 32| (0, 2, 12)(0, 1, 22)(0, 4, 24)(0, 16, 24) |
 | (0,  3, 22) | 32| (0, 6, 12)(0, 8, 16)(0, 3, 22)(0, 12, 24) |
 | (0,  5, 22) | 32| (0, 10, 12)(0, 5, 22)(0, 16, 24)(0, 20, 24) |
 | (0,  7, 22) | 32| (0, 12, 14)(0, 8, 16)(0, 7, 22)(0, 24, 28) |
 | (0,  9, 22) | 32| (0, 12, 18)(0, 9, 22)(0, 4, 24)(0, 16, 24) |
 | (0, 11, 22) | 32| (0, 8, 16)(0, 11, 22)(0, 12, 22)(0, 12, 24) |
 | (0, 13, 22) | 32| (0, 13, 22)(0, 16, 24)(0, 20, 24)(0, 12, 26) |
 | (0, 15, 22) | 32| (0, 8, 16)(0, 15, 22)(0, 24, 28)(0, 12, 30) |
 | (0, 17, 22) | 32| (0, 2, 12)(0, 17, 22)(0, 4, 24)(0, 16, 24) |
 | (0, 19, 22) | 32| (0, 6, 12)(0, 8, 16)(0, 19, 22)(0, 12, 24) |
 | (0, 21, 22) | 32| (0, 10, 12)(0, 21, 22)(0, 16, 24)(0, 20, 24) |
 | (0,  1, 23) | 16| (0, 2, 14)(0, 1, 23)(0, 8, 24)(0, 4, 28) |
 | (0,  2, 23) | 32| (0, 4, 14)(0, 8, 16)(0, 2, 23)(0, 8, 28) |
 | (0,  5, 23) | 16| (0, 10, 14)(0, 5, 23)(0, 8, 24)(0, 20, 28) |
 | (0,  6, 23) | 32| (0, 12, 14)(0, 8, 16)(0, 6, 23)(0, 24, 28) |
 | (0,  9, 23) | 16| (0, 14, 18)(0, 9, 23)(0, 8, 24)(0, 4, 28) |
 | (0, 10, 23) | 32| (0, 8, 16)(0, 14, 20)(0, 10, 23)(0, 8, 28) |
 | (0, 13, 23) | 16| (0, 13, 23)(0, 8, 24)(0, 14, 26)(0, 20, 28) |
 | (0, 14, 23) | 32| (0, 8, 16)(0, 14, 23)(0, 14, 28)(0, 24, 28) |
 | (0, 17, 23) | 16| (0, 2, 14)(0, 17, 23)(0, 8, 24)(0, 4, 28) |
 | (0, 18, 23) | 32| (0, 4, 14)(0, 8, 16)(0, 18, 23)(0, 8, 28) |
 | (0, 21, 23) | 16| (0, 10, 14)(0, 21, 23)(0, 8, 24)(0, 20, 28) |
 | (0, 22, 23) | 32| (0, 12, 14)(0, 8, 16)(0, 22, 23)(0, 24, 28) |
 | (0,  2, 25) | 32| (0, 4, 8)(0, 4, 18)(0, 16, 24)(0, 2, 25) |
 | (0,  3, 25) | 16| (0, 4, 12)(0, 6, 18)(0, 8, 24)(0, 3, 25) |
 | (0,  6, 25) | 32| (0, 12, 18)(0, 4, 24)(0, 16, 24)(0, 6, 25) |
 | (0,  7, 25) | 16| (0, 14, 18)(0, 8, 24)(0, 7, 25)(0, 4, 28) |
 | (0, 10, 25) | 32| (0, 4, 8)(0, 18, 20)(0, 16, 24)(0, 10, 25) |
 | (0, 11, 25) | 16| (0, 4, 12)(0, 18, 22)(0, 8, 24)(0, 11, 25) |
 | (0, 14, 25) | 32| (0, 4, 24)(0, 16, 24)(0, 14, 25)(0, 18, 28) |
 | (0, 15, 25) | 16| (0, 8, 24)(0, 15, 25)(0, 4, 28)(0, 18, 30) |
 | (0, 18, 25) | 32| (0, 4, 8)(0, 4, 18)(0, 16, 24)(0, 18, 25) |
 | (0, 19, 25) | 16| (0, 4, 12)(0, 6, 18)(0, 8, 24)(0, 19, 25) |
 | (0, 22, 25) | 32| (0, 12, 18)(0, 4, 24)(0, 16, 24)(0, 22, 25) |
 | (0, 23, 25) | 16| (0, 14, 18)(0, 8, 24)(0, 23, 25)(0, 4, 28) |
 | (0,  1, 26) | 32| (0, 4, 8)(0, 2, 20)(0, 16, 24)(0, 1, 26) |
 | (0,  3, 26) | 32| (0, 8, 12)(0, 8, 16)(0, 6, 20)(0, 3, 26) |
 | (0,  5, 26) | 32| (0, 8, 20)(0, 10, 20)(0, 16, 24)(0, 5, 26) |
 | (0,  7, 26) | 32| (0, 8, 16)(0, 14, 20)(0, 7, 26)(0, 8, 28) |
 | (0,  9, 26) | 32| (0, 4, 8)(0, 18, 20)(0, 16, 24)(0, 9, 26) |
 | (0, 11, 26) | 32| (0, 8, 12)(0, 8, 16)(0, 20, 22)(0, 11, 26) |
 | (0, 13, 26) | 32| (0, 8, 20)(0, 16, 24)(0, 13, 26)(0, 20, 26) |
 | (0, 15, 26) | 32| (0, 8, 16)(0, 15, 26)(0, 8, 28)(0, 20, 30) |
 | (0, 17, 26) | 32| (0, 4, 8)(0, 2, 20)(0, 16, 24)(0, 17, 26) |
 | (0, 19, 26) | 32| (0, 8, 12)(0, 8, 16)(0, 6, 20)(0, 19, 26) |
 | (0, 21, 26) | 32| (0, 8, 20)(0, 10, 20)(0, 16, 24)(0, 21, 26) |
 | (0, 23, 26) | 32| (0, 8, 16)(0, 14, 20)(0, 23, 26)(0, 8, 28) |
 | (0, 25, 26) | 32| (0, 4, 8)(0, 18, 20)(0, 16, 24)(0, 25, 26) |
 | (0,  1, 27) | 16| (0, 4, 12)(0, 2, 22)(0, 8, 24)(0, 1, 27) |
 | (0,  2, 27) | 32| (0, 8, 12)(0, 8, 16)(0, 4, 22)(0, 2, 27) |
 | (0,  5, 27) | 16| (0, 12, 20)(0, 10, 22)(0, 8, 24)(0, 5, 27) |
 | (0,  6, 27) | 32| (0, 8, 16)(0, 12, 22)(0, 12, 24)(0, 6, 27) |
 | (0,  9, 27) |  16| (0, 4, 12)(0, 18, 22)(0, 8, 24)(0, 9, 27) |
 | (0, 10, 27) | 32| (0, 8, 12)(0, 8, 16)(0, 20, 22)(0, 10, 27) |
 | (0, 13, 27) | 16| (0, 12, 20)(0, 8, 24)(0, 22, 26)(0, 13, 27) |
 | (0, 14, 27) | 32| (0, 8, 16)(0, 12, 24)(0, 14, 27)(0, 22, 28) |
 | (0, 17, 27) | 16| (0, 4, 12)(0, 2, 22)(0, 8, 24)(0, 17, 27) |
 | (0, 18, 27) | 32| (0, 8, 12)(0, 8, 16)(0, 4, 22)(0, 18, 27) |
 | (0, 21, 27) | 16| (0, 12, 20)(0, 10, 22)(0, 8, 24)(0, 21, 27) |
 | (0, 22, 27) | 32| (0, 8, 16)(0, 12, 22)(0, 12, 24)(0, 22, 27) |
 | (0, 25, 27) | 16| (0, 4, 12)(0, 18, 22)(0, 8, 24)(0, 25, 27) |
 | (0, 26, 27) | 32| (0, 8, 12)(0, 8, 16)(0, 20, 22)(0, 26, 27) |
 | (0,  2, 29) | 32| (0, 8, 20)(0, 16, 24)(0, 4, 26)(0, 2, 29) |
 | (0,  3, 29) | 16| (0, 12, 20)(0, 8, 24)(0, 6, 26)(0, 3, 29) |
 | (0,  6, 29) | 32| (0, 16, 24)(0, 20, 24)(0, 12, 26)(0, 6, 29) |
 | (0,  7, 29) | 16| (0, 8, 24)(0, 14, 26)(0, 20, 28)(0, 7, 29) |
 | (0, 10, 29) | 32| (0, 8, 20)(0, 16, 24)(0, 20, 26)(0, 10, 29) |
 | (0, 11, 29) | 16| (0, 12, 20)(0, 8, 24)(0, 22, 26)(0, 11, 29) |
 | (0, 14, 29) | 32| (0, 16, 24)(0, 20, 24)(0, 26, 28)(0, 14, 29) |
 | (0, 15, 29) | 16| (0, 8, 24)(0, 20, 28)(0, 15, 29)(0, 26, 30) |
 | (0, 18, 29) | 32| (0, 8, 20)(0, 16, 24)(0, 4, 26)(0, 18, 29) |
 | (0, 19, 29) | 16| (0, 12, 20)(0, 8, 24)(0, 6, 26)(0, 19, 29) |
 | (0, 22, 29) | 32| (0, 16, 24)(0, 20, 24)(0, 12, 26)(0, 22, 29) |
 | (0, 23, 29) | 16| (0, 8, 24)(0, 14, 26)(0, 20, 28)(0, 23, 29) |
 | (0, 26, 29) | 32| (0, 8, 20)(0, 16, 24)(0, 20, 26)(0, 26, 29) |
 | (0, 27, 29) | 16| (0, 12, 20)(0, 8, 24)(0, 22, 26)(0, 27, 29) |
 | (0,  1, 30) | 32| (0, 4, 24)(0, 16, 24)(0, 2, 28)(0, 1, 30) |
 | (0,  3, 30) | 32| (0, 8, 16)(0, 12, 24)(0, 6, 28)(0, 3, 30) |
 | (0,  5, 30) | 32| (0, 16, 24)(0, 20, 24)(0, 10, 28)(0, 5, 30) |
 | (0,  7, 30) | 32| (0, 8, 16)(0, 14, 28)(0, 24, 28)(0, 7, 30) |
 | (0,  9, 30) | 32| (0, 4, 24)(0, 16, 24)(0, 18, 28)(0, 9, 30) |
 | (0, 11, 30) | 32| (0, 8, 16)(0, 12, 24)(0, 22, 28)(0, 11, 30) |
 | (0, 13, 30) | 32| (0, 16, 24)(0, 20, 24)(0, 26, 28)(0, 13, 30) |
 | (0, 15, 30) | 32| (0, 8, 16)(0, 24, 28)(0, 15, 30)(0, 28, 30) |
 | (0, 17, 30) | 32| (0, 4, 24)(0, 16, 24)(0, 2, 28)(0, 17, 30) |
 | (0, 19, 30) | 32| (0, 8, 16)(0, 12, 24)(0, 6, 28)(0, 19, 30) |
 | (0, 21, 30) | 32| (0, 16, 24)(0, 20, 24)(0, 10, 28)(0, 21, 30) |
 | (0, 23, 30) | 32| (0, 8, 16)(0, 14, 28)(0, 24, 28)(0, 23, 30) |
 | (0, 25, 30) | 32| (0, 4, 24)(0, 16, 24)(0, 18, 28)(0, 25, 30) |
 | (0, 27, 30) | 32| (0, 8, 16)(0, 12, 24)(0, 22, 28)(0, 27, 30) |
 | (0, 29, 30) | 32| (0, 16, 24)(0, 20, 24)(0, 26, 28)(0, 29, 30) |
 | (0,  1, 31) | 16| (0, 8, 24)(0, 4, 28)(0, 2, 30)(0, 1, 31) |
 | (0,  2, 31) | 32| (0, 8, 16)(0, 8, 28)(0, 4, 30)(0, 2, 31) |
 | (0,  5, 31) | 16| (0, 8, 24)(0, 20, 28)(0, 10, 30)(0, 5, 31) |
 | (0,  6, 31) | 32| (0, 8, 16)(0, 24, 28)(0, 12, 30)(0, 6, 31) |
 | (0,  9, 31) | 16| (0, 8, 24)(0, 4, 28)(0, 18, 30)(0, 9, 31) |
 | (0, 10, 31) | 32| (0, 8, 16)(0, 8, 28)(0, 20, 30)(0, 10, 31) |
 | (0, 13, 31) | 16| (0, 8, 24)(0, 20, 28)(0, 26, 30)(0, 13, 31) |
 | (0, 14, 31) | 32| (0, 8, 16)(0, 24, 28)(0, 28, 30)(0, 14, 31) |
 | (0, 17, 31) | 16| (0, 8, 24)(0, 4, 28)(0, 2, 30)(0, 17, 31) |
 | (0, 18, 31) | 32| (0, 8, 16)(0, 8, 28)(0, 4, 30)(0, 18, 31) |
 | (0, 21, 31) | 16| (0, 8, 24)(0, 20, 28)(0, 10, 30)(0, 21, 31) |
 | (0, 22, 31) | 32| (0, 8, 16)(0, 24, 28)(0, 12, 30)(0, 22, 31) |
 | (0, 25, 31) | 16| (0, 8, 24)(0, 4, 28)(0, 18, 30)(0, 25, 31) |
 | (0, 26, 31) | 32| (0, 8, 16)(0, 8, 28)(0, 20, 30)(0, 26, 31) |
 | (0, 29, 31) | 16| (0, 8, 24)(0, 20, 28)(0, 26, 30)(0, 29, 31) |
 | (0, 30, 31) | 32| (0, 8, 16)(0, 24, 28)(0, 28, 30)(0, 30, 31) |


<br>

------

References and Footnotes
------

[^integers]:   Mentally replace that by $\mathbb{Z}_{2^b}$ or $\mathbb{Z}/2^{b}$ if you're offended.
[^modinverse]: **"Integer multiplicative inverse via Newton's method"**, 2017 ([page]({{site.base}}/math/2017/09/18/ModInverse.html))
[^pscarab]:    **"Hash functions"**, Bret Mulvey ([page](http://papa.bretmulvey.com/post/124027987928/hash-functions))
[^pixar]:      **"Correlated Multi-Jittered Sampling"**, Andrew Kensler 2013 ([PDF](http://graphics.pixar.com/library/MultiJitteredSampling/paper.pdf))
[^permmatrix]: Wikipedia: Permutation matrix ([page](http://en.wikipedia.org/wiki/Permutation_matrix))
[^simp]:       Well $\mathbb{F}_{2}$ is more general than math on bits...don't need any of that here.
[^shiftm]:     Wikipedia: Shift Matrix ([page](http://en.wikipedia.org/wiki/Shift_matrix))
[^bijection]:  Wikipedia: Bijection ([page](http://en.wikipedia.org/wiki/Bijection))
[^exchange]:   Wikipedia: Exchange matrix ([page](http://en.wikipedia.org/wiki/Exchange_matrix))
[^cayley]:     Wikipedia: Cayley's theorem ([page](http://en.wikipedia.org/wiki/Cayley%27s_theorem))
[^qsymbol]:    Wikipedia: Q-Pochhammer ([page](http://en.wikipedia.org/wiki/Q-Pochhammer_symbol))
[^arndt]:      **"Matters Computational"**, Jorg Arndt, 2010 ([page](http://www.jjj.de/fxt/))
[^circulant]:  Wikipedia: Circulant matrix ([page](https://en.wikipedia.org/wiki/Circulant_matrix#Properties))
[^composition]:  A composition of bijections is a bijection.
[^invmix]:      **"Can an involution be a competitive bit finalizer?"**, 2019 ([page]({{site.base}}/math/2019/08/20/InvFinalizer.html))


<script>

// can't be bothered
function popcount(x) {var m=1,r=0; while (m>0) {if(x&m){r++;} m = (m<<1)&0xffffffff;} return r;}

function grid_draw_d(ctx,x,y,len)
{
  for (var i=0; i<len; i++) { ctx.fillRect(x, y, 4, 4);	x += 5;	y += 5; }
}

function grid_draw_ld(ctx,pos) { grid_draw_d(ctx,1,5*pos+1,32-pos); }
function grid_draw_ud(ctx,pos) { grid_draw_d(ctx,5*pos+1,1,32-pos); }

function grid_draw_xorrot(grid,k)
{
  k = k & 0xffffffff;
  var ctx = grid.ctx;
  var pos = 1;
  if ((k & 1)==1) { grid_draw_ld(ctx,0); } k = k >>> 1;
  while(k != 0) {
    if ((k & 1)==1) { grid_draw_ld(ctx, pos); grid_draw_ud(ctx, 32-pos); }
	pos++; k = k >>> 1;
  }
}



function grid_draw(grid)
{
  var ctx  = grid.ctx;
  var size = 5*grid.dim+1;

  ctx.fillStyle = grid.bgfill;
  ctx.fillRect(0, 0, size, size);
  ctx.fillStyle = grid.fgfill;
  grid_draw_xorrot(grid, grid.k);
  grid.cnt++;
  
  if (grid.cnt == 10) {
    var p;
	grid.k = (grid.k + 0x3504f333) & 0xFFFFFFFF;
	
    do { grid.k++; p = popcount(grid.k); } while(p < 3 ||  (p&1)!=1);

	var text = "(";
	var k = grid.k;
	var p = 0;
	var x = false;
	
	do { 
	  if (k & 1) { if (x) {text=text+","; } text = text + p; x=true;}
	  p++; k = k >>> 1;
	} while(k != 0);


	document.getElementById("dfig").textContent = text + ")";
    grid.cnt=0; 
  }
  
  window.setTimeout(grid.draw, 20);
}

function grid_make(name)
{
  var c = document.getElementById(name);
  var x = c.getContext('2d');

  var grid = {
	ctx:    x,
	draw:   null,
	k:      3,
	dim:    32,
	cnt:    0,
	bgfill: 'rgba(230,240,255, 0.3)',
	fgfill: '#f00'
  };

  grid.draw = function() { grid_draw(grid); };

  return grid;
}

var grid = grid_make('fig');

grid.draw();

</script>
