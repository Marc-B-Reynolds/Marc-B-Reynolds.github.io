---
layout:       post
title:        'FMA: solve quadratic equation'
categories:   [math]
tags:         [FMA,floating point]
description:  Micro-post on title.
plotly:       false
---

\\
Micro-post: Solve the [quadratic equation](https://en.wikipedia.org/wiki/Quadratic_equation) without [catastrophic cancellation](https://en.wikipedia.org/wiki/Loss_of_significance) using FMA and no branching.

Given $a$, $b$ and $c$ we want to find one or both roots $r_0$ and $r_1$:

$$ ax^2 + bx + c = a\left(\left(x-r_0\right)\left(x-r_1\right)\right) $$

\\
Recall we have two sources of catastrophic cancellation:

1. $b^2 \gg \left\\| 4ac \right\\| $ which we'll fix with an algebraic rewrite
2. $b^2 \approx 4ac $   which we'll fix by computing in extended precision.

Calling $r_1$ the larger and $r_0$ the smaller magnitude root we can rewrite the standard equation as follows:

$$ \begin{align*}
r_1 & = \frac{-\left(b + \text{sgn}\left(b\right)\sqrt{b^2-4ac}\right) }{2a} \\
    & = \frac{b + \text{sgn}\left(b\right)\sqrt{b^2-4ac}}{-2a} \\
\\
r_0 & = \frac{2c}{-\left(b + \text{sgn}\left(b\right)\sqrt{b^2-4ac}\right) } \\
    & = \frac{-2c}{b + \text{sgn}\left(b\right)\sqrt{b^2-4ac}}
\end{align*} $$

\\
where $\text{sgn}$ is the sign function:

$$
\text{sgn}\left(x\right) =
\begin{cases}
1  & x \geq 0 \\[2ex]
-1 & x < 0
\end{cases}
$$

\\
The remaining part is computing $b^2-4ac$ in extended precision which is a special case of extended precision $ab+cd$. This is covered in detail in [*"Further analysis of Kahan's algorithm for the accurate computation of 2x2 determinants"*](https://hal.inria.fr/ensl-00649347/en)


\\
Toy code which assumes two real roots (including $r_0=r_1$)

{% highlight c %}

// Please do use -fno-math-errno
inline float f32_sqrt(float x) { return sqrtf(x); }

typedef struct { float h,l; } f32_pair_t;

// ab-cd
// * r within +/- 3/2 ulp
// * r within +/-   1 ulp, if sign(ab) != sign(cd)
inline float f32_mms(float a, float b, float c, float d)
{
  float t = c*d;
  float e = fmaf(c,d,-t);
  float f = fmaf(a,b,-t);
  return f-e;
}

// larger magnitude root
float f32_quadratic_max(float a, float b, float c)
{
  float t0 = f32_sqrt(f32_mms(b,b,4.f*a,c));
  float t1 = b+copysignf(t0,b);
  return t1/(-2.f*a);
}

// smaller magnitude root
float f32_quadratic_min(float a, float b, float c)
{
  float t0 = f32_sqrt(f32_mms(b,b,4.f*a,c));
  float t1 = b+copysignf(t0,b);
  return (-2.f*c)/t1;
}

void f32_quadratic(f32_pair_t* r, float a, float b, float c)
{
  float t0 = f32_sqrt(f32_mms(b,b,4.f*a,c));
  float t1 = b+copysignf(t0,b);

  r->h = t1/(-2.f*a);
  r->l = (-2.f*c)/t1;
}

{% endhighlight %}

