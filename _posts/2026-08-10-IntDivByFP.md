---
layout:       post
title:        Moving integer division to floating-point is trivial
categories:   [math]
tags:         [integer]
description:  "Shows that floating point trunc(x/y) equals integer x/y provided the integers fit in the floating point format."
plotly:       false
---

\\
Integer division `q=(x/y)` and remainder (of Euclidean division) `r=(x%y)` hardware operations are very sad
on current hardware. Typically very long latency and poor throughput. In contrast floating-point division
is pretty happy: shorter latency, higher throughput and often more execution units to perform the
operation. So there are cases it could be interesting to move some integer div/mod operations to floating point.
But it's PITA right?  Actually I think it's easy. The math is pretty straightforward so if I've made a mistake
I expect to find out rather soon.

My claim is: for two integers `x` & `y` (signed or unsigned) that fit in 53/24 bits for double/single
precision respectively, with both promoted to floating point then:

```c?encoding=ascii
// floating-point:
//   d is the same integer as integer divide    x / y
//   m is the same integer as integer remainder x % y

d = trunc(x/y);     // floor works for unsigned
m = -fma(d,y,-x);   // fma required

// NOTE: if only want 'd' and it's being converted to an
// integer then the truncate or floor operation is
// free in the float to integer conversion.
```

<br>
in standard rounding mode. 

From here on I will only consider unsigned integers since it's the harder case (signed have a smaller max magnitude) but recall
that floating point effective stores a signed magnitude quantity.

Some practical points:
* this works by setting rounding-mode to `TOWARD_ZERO` and then back once done but hitting control words is often very expensive
* some hardware have division opcodes that allow choosing the rounding mode on some operations like division. 
  Example *AVX-512* has instrinsic: `_mm_div_round_sd` with example listed latency of 14 cycles.
* converting int to float and back has a cost. (NOTE: the majority of x64 CPUs don't have an op for unsigned to float and
  back. So this is an implemenation concern.)
* obviously the floating point precision $p$ places the upper bound on integer width. If this works for $p$ then all smaller
  widths trivally work. Examples 32-bit integers in doubles ($p=53$) and 16-bit integers in singles ($p=24$). In this case it's
  possible that there's some "trickery" to side-step the int to float.
* obviously the most promising case is for working in SIMD to amortize overhead.

Note that there are proven methods (SEE: *[Formally verified 32- and 64-bit integer division using double-precision floating-point arithmetic](https://arxiv.org/abs/2207.08420)*)

The standard rounding mode is: *round-to-nearest (ties to even)* and the first important observation is the *ties*
part. A tie happens when the exact result of an operation is exactly at the *midpoint* between two 
floating point numbers. Floating point division (in our case of base-2 and same working precision) has zero
midpoints (SEE: *[Midpoints and exact points of some algebraic functions in floating-point arithmetic](https://ens-lyon.hal.science/ensl-00409366/)* section 6.1, corollary 1).
Therefore there's never a *tie* and only *round-to-nearest* portion applies.

Let's look an example using a 4-bit precision floating point format with all possible configurations where the exact 
result is less than one but as close as possible to rounding up:
      
          |GRS            x = don't care
     .1111|0xx   .1111    (no rounding)
     .1111|100            tie case is impossible
     .1111|101  1.000     (round up)
     .1111|11x  1.000     (round up)

<br>
Prior to rounding the hardware computes three extra digits: *guard bit* (G), *round bit* (R) and *sticky bit* (S) but since
the *tie* case is impossible we only need to know *G*.  Therefore in any format where the exact result has $r$ bits for the
fractional part then for rounding to the next integer to occur requires the fractional part to have at least $f+1$ leading ones.
I'm claiming that this is impossible with legal inputs.

\\
Since we only need to consider what happens with the fractional part let's breakdown the exact result of $x/y$ into its integer $n$ and fractional parts:

$$
\frac{x}{y} = n + \frac{a}{b}
$$

\\
(where $b=y$ just to be less of an eyesore) so obviously:

$$
\frac{a}{b} \in \left[0,1\right)
$$


Given a $p$ precision binary floating point format the division produces a $d$-bit integer leaving $r$-bits for the remainder:

$$
\begin{align*}
d & =  \left \lfloor \log_2 \left(n \right)+1 \right \rfloor \\
r & = p-d    \\
\end{align*}
$$

For a $r$-bit fractional part we can define the lower and upper bound on it's range of value for $b$:

$$
\begin{align*}
\func{b_l}{r} & = 2^r \\
\func{b_u}{r} & = 2^{r+1}-1
\end{align*}
$$

At this point we're pretty much done with $p$, $n$ and its derived $d$. We want to know the maximum fractional part possible for $r$-bits.

The closest value to $1$ we can create with the lower is:

$$
\begin{align*}
\frac{2^r-1}{2^r} = \sum_{i=1}^r \frac{1}{2^i}
\end{align*}
$$

\\
Rewritten as the sum tells us that this is the smallest number that produces $r$ leading ones in it's binary expansion.

Now repeating with the upper bound's closest to $1$ value and compare to the smallest number with $r+1$ leading ones (the lower bound of $r+1$):

$$
\begin{align*}
\frac{2^{r+1}-2}{2^{r+1}-1} &< \frac{2^{r+1}-1}{2^{r+1}-0}
\end{align*}
$$

\\
where the minus zero is to align the expressions and highlight the pattern. Anyway this inequality means that the upper bound is limited to having $r$ leading ones.

We can actually be specific about the binary expansion of bounds. The lower is $r$ ones and terminating. The upper's exact value is purely cyclic with a period of $r+1$ with $r$ ones followed by one zero.  Here's a table for the first four:

{: .center }
|$r$  | lower         | binary expansion |upper          | binary expansion |
|:---:|   :---:       | :---             |  :---:        | :---             |
| 1   |$\frac{1}{2}$  |  `.1`            |$\frac{2}{3}  $|  `.(10)`         |
| 2   |$\frac{3}{4}$  |  `.11`           |$\frac{6}{7}  $|  `.(110)`        |
| 3   |$\frac{7}{8}$  |  `.111`          |$\frac{14}{15}$|  `.(1110)`       |
| 4   |$\frac{15}{16}$|  `.1111`         |$\frac{30}{31}$|  `.(11110)`      |


\\
But anyway the important part:

<div class="alert alert-success" role="alert" markdown="1">
For legal input: $\frac{x}{y}$ cannot produce a fractional part large enough to cause rounding to the
next integer so the integer part of the floating point result is identical to integer result.
</div>

<br>
Let's briefly tie the parts together by considering some divisors.

Dividing by:
* $1$ is correct
* $2$ only modifies the exponent and no rounding of any kind occurs.
* $3$ which is $r=1$, the largest `x % 3` value is $2$ so max fraction is $\frac{2}{3}$. Can't round to next integer.
* $\set{4,5,6,7}$ have $r=2$. The exact divisor doesn't matter because they all start with `.110x`. 

And so on all possible divisors.



