---
layout:       post
title:        "Special cased quaternion log and exp"
categories:   [quaternions]
tags:         [animation, rotation]
description:  Describes efficient implementations of scaled quaternion log and exp.
---

\\
I'm going to quickly run through how to implement efficient and accurate specialized versions of `log`
and `exp` of unit quaternions where the exponential space is folded to the unit ball. The
goal is to have an error bound no greater than a straight forward implementation using standard
functions. Although I'm only showing scalar these are designed to be *SIMD* friendly.  The specialized
functions are:

$$
\begin{align*}
\mathbf{v} & = \frac{2}{\pi}\tfunc{log}{\abs{Q}}    \\
Q          & = \tfunc{exp}{\frac{\pi}{2}~\mathbf{v}}
\end{align*}
$$ 

\\
where $\abs{Q}$ is shorthand for canonicalize to positive scalar. (negate $Q$ if $w$ is negative). The
impact of this specialization is that the exponential space is limited to torque minimal angles.


<br>

------

Linear to exponential space
------

\\
First let's look at an accurate version of standard `log` using standard library functions.

Given an arbitrary quaternion $Q$ with implied scale factor $m$ and unit bivector $\mathbf{u}$:

$$
\begin{align*}
Q & = m~\cos\left(\theta\right)+m~\sin\left(\theta\right)~\mathbf{u} \\
  & = w + \mathbf{v}
\end{align*}
$$ 

\\
then we want to compute:

$$
\tfunc{log}{Q} = \theta~\mathbf{u}
$$ 


```c?encoding=utf8
// log(q) : standard functions 
static inline vec3f_t quatf_log(quatf_t q)
{
  float x = q[3];               // m cos(Θ)
  float t = quat_bnorm_fma(q);  // m² sin(Θ)²
  float y = sqrtf(t);           // |V| = m sin(Θ)
  float a = atan2f(y,x);        // Θ : this does all the heavy lifting
  float s = a/y;                // scale factor

  // 'y' can only be "normal" or zero. In the limit the scale factor 's'
  // approaches one (as y approaches zero). Alternately this could set
  // s to zero for y==0 which would filter tiny angles to the identity.
  s = (y != 0.f) ? s : 1.f;

  return s*quat_bivector(q);
}
```

<br>
and simply wrapping the previous function to get our specialized version and note
the approximate distance error (which is roughly a half turn angle error).

```c?encoding=utf8
// wrap log:  distance error ~2.54740144e-07
static inline vec3f_t quatf_fem_hq(quatf_t q)
{
  static const float K = 0x1.45f306p-1f;  // 2/π

  q = (q[3] >= 0.f) ? q : -q;

  return K*quat_log(q);
}
```

<br>
The goal isn't to be error bound competitive with the previous. Instead we'll build
version (using standard functions) that assume $Q$ is exactly unit magnitude. I'm going
to skip on showing a standard `log` version and jump directly to the specialized:

```c?encoding=utf8
// (2/π) log(|q|) : standard functions: move to acos and computing |V| as sqrt(1-w²)
static inline vec3f_t quatf_fem_std(quatf_t q)
{
  static const float K = 0x1.45f306p-1f;  // 2/π

  float k = copysignf(K,q[3]);            // sgn(q.w) 2/π
  float w = fabsf(q[3]);                  // |q.w|

  w = (w <= 1.f) ? w : 1.f;               // clamp |w| to range
                                          // distance error (approximate)
  float y = sqrtf(-fmaf(w,w,-1.f));       //   3.40125329e-07
//float y = sqrtf((1.f-w)*(1.f+w));       //   3.63783706e-07 ('y' error bound is same as FMA above)
//float y = sqrtf(1.f-w*w);               //   5.11455530e-07
//float y = sqrtf(quat_bnorm_fma(q));     //   4.68131871e-06
  float a = acosf(w);
  float s = a/y;
  
  s  = (y != 0.f) ? s : 1.f;              // same comments on limit
  s *= k;

  return s*quat_bivector(q);
}
```

<br>
Notice I have a series of commented-out computations for $y = \abs{\mathbf{v}}$. Throughout this post I'm assuming all targets have hardware `FMA`.
If we need to hit a target without then for all other cases we can simple replace with a multiply and add (note for bit indentical results across targets
then this choice needs to be consistent). Elsewhere this will only result in a small increase in error bound. For the computation of $y$ then we want
to use the first commented out version.  The next is to show that swapping out the *FMA* into a multiply/add sequence for $y$ is very sad. And the last
is to show that the accurate way to compute $y$ ends up being the worst choice since that information is not being used to compute the angle. If we
modified `quatf_log` to instead compute `float y = sqrtf(-fmaf(w,w,-1.f))` then we run into the same problem and the error our wrapped `quatf_fem_hq` would
increase to `5.11220598e-07` making it worse than the above version.

Notice that to create a direct implementation we simply need to formulate an approximation for the scaling factor $s$ above which is solely a function of $w$:

$$
\func{s}{w} = \frac{2}{\pi}~\frac{\tfunc{acos}{w}}{\sqrt{1-w^2}}
$$

\\
There are a number of ways we approximate this and the [file](https://github.com/Marc-B-Reynolds/Stand-alone-junk/blob/main/src/SFH/f32_math/internal/quatf_fem.h)
(currently) has 19 choices covering a range of error bounds using three different forms: a straight polynomial, a rational and a transformed approximation.  The latter ends up being
uninteresting without an accurate hardware reciprocal square root (IMHO).  Here I'm going to show the first rational approximation that is has a smaller error bound than
the above.

```c?encoding=utf8
// ~3.37359953e-07
static inline vec3f_t quatf_fem(quatf_t q)
{
  static const float C   = 0x9.81b88p-4f;
  static const float P[] = {0x9.f2a54p-4f, 0x2.2b9cfp-4f, 0x7.32ff48p-12f};
  static const float Q[] = {0x7.9a5a8p-4f, 0xe.7215p-8f};
  
  float w = fabsf(q[3]);
  float n = P[2];
  float d = Q[1];

  n = fmaf(n,w,P[1]);
  n = fmaf(n,w,P[0]); 
  n = fmaf(n,w,C);
  n = copysignf(n,q[3]);
  
  d = fmaf(d,w,Q[0]);
  d = fmaf(d,w,1.f);
  d = fmaf(d,w,C);
  
  return (n/d)*quat_bivector(q);
}
```

<br>
As noted in the linked file the above rational approximation was constructed using [rminimax](https://gitlab.inria.fr/sfilip/rminimax) with the following command line:

<pre><small>
ratapprox --function="2*(acos(x)/(sqrt(1-x^2)))/pi" --dom=[0,0.99999999] --denF=[SG] --numF=[SG] --num=[1,x,x^2,x^3] --den=[1,x,x^2,x^3] --weight=1 --output=fem_33.sollya
</small></pre>



<br>

------

Exponential to linear space
------

\\
Our special case inverse mapping becomes:  given $\mathbf{v} = \theta~\mathbf{u}$

$$
\begin{align*}
Q & = \tfunc{exp}{\frac{\pi}{2}~\mathbf{v}} \\
  & = \tfunc{cos}{\frac{\pi}{2}~\theta} + \tfunc{sin}{\frac{\pi}{2}~\theta}~\mathbf{u} \\
  & = \tfunc{cos}{\frac{\pi}{2}~\theta} + \frac{\tfunc{sin}{\frac{\pi}{2}~\theta}}{\sqrt{\mathbf{v} \cdot \mathbf{v}}}~\mathbf{v}
\end{align*}
$$

<br>
So a straight-forward version with standard functions:

```c?encoding=utf8
// exp(π/2 V) special cased for 'V' in unit ball
//  given V = ΘU (in unit ball, U = unit bivector)
//  returns cos(π/2 Θ) + sin(π/2 Θ) U
//  distance error: ~2.63713986e-07
static inline quatf_t quatf_iem_naive(vec3f_t v)
{
  float d = vec3_norm_fma(v);        // v∙v
  float x = sqrtf(d);                // |v|
  float a = ((float)(0.5*M_PI))*x;   // this casted double happens to be correctly rounded in singles
  float w = cosf(a);                 // or sincospi if you have it
  float s = sinf(a)/x;
  
  s = (x != 0.f) ? s : 1.f;          // correct any limit case
  
  return quatf_bs(s*v,w);            //   since this is not for real use
}
```

<br>
To construct a direct version let's note that $\cos$ is even and can be approximated by:

$$
\tfunc{cos}{\frac{\pi}{2}~\theta} = 1 + \theta^2~\func{P}{\theta^2}
$$

\\
where $P$ is some polynomial and that computing $\theta$ starts as:

$$
d = \mathbf{v} \cdot \mathbf{v} = \theta^2 
$$

\\
Plugging that back into $\cos$, noting that $\theta$ is positive and let $p = \func{P}{d}$:

$$
\begin{align*}
\tfunc{cos}{\frac{\pi}{2}~\theta} & = 1 + dp \\
\tfunc{sin}{\frac{\pi}{2}~\theta} & = \sqrt{1-\left(1 + dp\right)^2} \\
                                  & = \sqrt{-2dp-d^2~p^2}
\end{align*}
$$

\\
which gives us:

$$
\begin{align*}
\frac{\tfunc{sin}{\frac{\pi}{2}~\theta}}{\sqrt{d}} & = \frac{\sqrt{-2dp-d^2~p^2}}{\sqrt{d}} \\
                                                   & = \sqrt{\frac{-2dp-d^2~p^2}{d}}        \\
                                                   & = \sqrt{-2p-d~p^2}
\end{align*}
$$

\\
Like for the forward function here's a [file](https://github.com/Marc-B-Reynolds/Stand-alone-junk/blob/main/src/SFH/f32_math/internal/quatf_iem.h) with a selection of implementations (currently 3) and I've choosen to show the minimum that outperforms the error bound of the above.

```c?encoding=utf8
// ~2.30486819e-07
static inline quatf_t quatf_iem(vec3f_t v)
{
  // cos(π/2 Θ) ≈ 1 + Θ²P(Θ²) where Θ² = d = v∙v
  //   coefficients of polynomial P
  static const float C[] = {0x1.c2a574p-11f, -0x1.5502d4p-6f,
                            0x1.03bd9ep-2f,  -0x1.3bd3bp0f};

  // computation of 's'
  // sin(π/2 Θ)/Θ = sqrt(1-w²)/Θ        : from cosine
  //              = sqrt((1-w²)/d)      : pull Θ inside sqrt
  //              = sqrt((1-(1+dp)²)/d) : w = 1+dp
  //              = sqrt(-dp²-2p)       : reduce
  float d = vec3_norm_fma(v);           // v∙v
  float p = f32_horner_3(d,C);          // p = P(d)
  float s = sqrtf(-fmaf(d,p*p,p+p));    // sin(π/2 Θ)/Θ  *
  float w = fmaf(d,p,1.f);              // cos(π/2 Θ)    *

  // { sin(π/2 Θ)/Θ v, cos(π/2 Θ) }
  return quatf_bs(s*v,w);
}
```

<br>
The two lines with asterisks are using FMAs structured to lower the error (which is why $pd$ is being computed twice) and I haven't bother to think through the best way to compute these quantities without FMA.


<br>

------

Success?
------

\\
Let's look at how successful the choosen example routines are performing by generating uniform random unit quaternions,
map to exponential space, back to linear and measuring the 3D rotation angle verses the original:


```c?encoding=utf8
// computes the 3D rotation (instead of intrinsic) angle between
// A & B measured in radians.
static inline double quatf_relative_angle_hq(quatf_t A, quatf_t B)
{
  quatd_t a = quatf_promote(A);               // promote to doubles
  quatd_t b = quatf_promote(B);               //   ditto
  quatd_t r = quatd_mul(a,quat_conj(b));      // compute relative: AB^*

  // use atan2 for higher accuracy
  return 2.0*atan2(sqrt(vec3_dot_fma(r,r)),r[3]);
}
```

<br>

and we'll repeat this with a set of forward/inverse functions:

1. reference versions (which are promote to double and standard function. see linked files)
2. perturb. No mappings: simply flips the lowest bit of each of the four components.
3. standard function versions
4. the choosen direct implementations contained in this post


|         fem   |            iem| radians       | example |
|  :---:        | :---:         |  :---:        |:---:    |
| quatf_fem_ref | quatf_iem_ref |2.07116297e-07 | `{ 0x1.1e561cp-1,-0x1.259016p-1,-0x1.2a5062p-1, 0x1.1a6be0p-3}` |
| perturb       | N/A           |2.13446408e-07 | `{ 0x1.3ea590p-2, 0x1.450110p-1, 0x1.0011cep-1, 0x1.000c62p-1}` |
| quatf_fem_std | quatf_iem_std |1.06645021e-06 | `{-0x1.27ba00p-5, 0x1.5eec0ap-1,-0x1.73f2a0p-1, 0x1.1a20c0p-5}` |
| quatf_fem     | quatf_iem     |1.07541018e-06 | `{-0x1.1f4bd0p-3, 0x1.104f90p-1, 0x1.ab9594p-1, 0x1.eeb800p-9}` |

\\
So we have:
* the reference pair introduces less error than simply perburbing
* our example pair is slightly worse than the standard function versions (abs difference is less that $9 \times 10^{-9}$) but good enough IMHO to call a success


