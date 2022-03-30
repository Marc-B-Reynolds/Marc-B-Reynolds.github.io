---
layout: post
title:      "Quaternion swing twist decomposition"
categories: [quaternions]
tagline:    "factor with respect to a given axis and the plane orthogonal to it"
tags :      [rotations, decomposition]
description : 'Given a unit quaternion $Q$ factor into a pair such that $Q = Q_rQ_a$ where $a$ is a predefined axis and $r$ is an orthogonal axis.  This is commonly called swing-twist decomposition.'
plotly: true
---

<div class="alert alert-success" role="alert" markdown="1">
<b>TL;DR</b> version: If you just want code snippets and enough info to code it up then just jump to *Comments on implementations*. 
</div>

This is expanding on a detail from a previous post ([*"Hopf coordinate conversion and torque minimal factorization"*](https://marc-b-reynolds.github.io/quaternions/2017/05/12/HopfCoordConvert.html)) which is that a quaternion $Q$ can be factored into:

$$ \begin{equation} \label{def}
Q = Q_sQ_t
\end{equation} $$

\\
which is commonly referred to as *"swing-twist decomposition"* where $Q_t$ (the twist) is a rotation about a pre-chosen axis and $Q_s$ (the swing) is some orthogonal axis. There is example code for doing this (`quat_to_tmtwist`) but it is storing the results in a specialized structure that contains both factors and skips on the zero components. The inverse function (`tmtwist_to_quat`) is just a specialized product with this knowledge.

So let's skim: given a unit quaternion:

$$ \begin{equation} \label{Qr}
Q = w + \left(x,~y,~z\right) 
\end{equation} $$ 

\\
we can factor into pair with respect to any rotation axis (the *twist*). For the three "coordinate" axes we have:

$$ 
\definecolor{foo}{rgb}{0.8, 0.2, 0.2}

\begin{align}
Q_{tx} &= \frac{w}{\sqrt{w^2+x^2}} + \left(\frac{x}{\sqrt{w^2+x^2}},~0,~0\right) 
& Q_{sx} &= \sqrt{w^2+x^2} + \left(0,~\frac{wy{\color{foo}-}xz}{\sqrt{w^2+x^2}},~\frac{wz{\color{foo}+}xy}{\sqrt{w^2+x^2}}\right) \label{fx} \\

Q_{ty} &= \frac{w}{\sqrt{w^2+y^2}} + \left(0,~\frac{y}{\sqrt{w^2+y^2},~0}\right) 
& Q_{sy} &= \sqrt{w^2+y^2} + \left(\frac{wx{\color{foo}+}yz}{\sqrt{w^2+y^2}},~0,~\frac{wz{\color{foo}-}xy}{\sqrt{w^2+y^2}}\right) \label{fy} \\

Q_{tz} &= \frac{w}{\sqrt{w^2+z^2}} + \left(0,~0,~\frac{z}{\sqrt{w^2+z^2}}\right) 
& Q_{sz} &= \sqrt{w^2+z^2} + \left(\frac{wx{\color{foo}-}yz}{\sqrt{w^2+z^2}},~\frac{wy{\color{foo}+}xz}{\sqrt{w^2+z^2}},~0\right) \label{fz} \\

\end{align}
$$


\\
Where $Q_{ta}$ is the *twist* about $\mathbf{a}$ and $Q_{sa}$ is the *swing* about an axis orthogonal to $a$. 

Now all of these factorizations have an ill conditioned region. When $t$ is approaching zero the input is approaching the maximum magnitude rotation ($\pi$) and axis is likewise approaching being independent of one one we're attempting to factor (again $\mathbf{z}$ for the working example).

The order of composition can also be reverse:

$$ Q = Q_tQ_r$$

The twist expressions remain the same and the swings are simply negate the component of the twist ($x$ is negate for twist about $\mathbf{x}$ and so forth) when the twist axis is coordinate frame aligned.

$$ 
\begin{align}
 Q_{rx} &= \sqrt{w^2+x^2} + \left(0,~\frac{wy{\color{foo}+}xz}{\sqrt{w^2+x^2}},~\frac{wz{\color{foo}+}xy}{\sqrt{w^2+x^2}}\right) \label{rx} \\
 Q_{ry} &= \sqrt{w^2+y^2} + \left(\frac{wx{\color{foo}-}yz}{\sqrt{w^2+y^2}},~0,~\frac{wz{\color{foo}+}xy}{\sqrt{w^2+y^2}}\right) \label{ry} \\
 Q_{rz} &= \sqrt{w^2+z^2} + \left(\frac{wx{\color{foo}+}yz}{\sqrt{w^2+z^2}},~\frac{wy{\color{foo}-}xz}{\sqrt{w^2+z^2}},~0\right) \label{rz} \\
\end{align}
$$

The derivations are given below. The factorization has the following properties (using $\mathbf{z}$ as an example): both $Q$ and $Q_{sz}$ rotate $\mathbf{z}$ (and the z component of a vector) to the same value:

$$ \mathbf{z}' = Q~\mathbf{sz}~Q^* = Q_{xy}~\mathbf{z}~Q_{sz}^*$$

and $Q_{tz}$ leaves $\mathbf{z}$ (and z components) unmodified:

$$  Q_{tz}~\mathbf{z}~Q_{tz}^* = \mathbf{z} $$

Both $Q$ and $-Q$ factorize to the same *swing*, the *twist* is negated and the reconstruction (in absence of filtering) produces the input.

Rewritting the factors using two temporary variables we have:

<div class="container">
  <div class="row">
    <div class="col-sm">
$$ 
\begin{array}{rl}
t &=  \sqrt{w^2+z^2} \\
s &=  1/t           \\
Q_{tz} &= sw + \left(0,~0,~sz\right) \\
Q_{sz} &= t  + s\left(\left(wx-yz\right),~\left(wy+xz\right),~0\right) \\

\end{array}
$$

    </div>
    <div class="col-sm">
$$ 
\begin{array}{rl}
t &=  w^2+z^2 \\
s &=  1/\sqrt{t} \\
Q_{tz} &= sw + \left(0,~0,~sz\right) \\
Q_{sz} &= st + s\left(\left(wx-yz\right),~\left(wy+xz\right),~0\right) \\

\end{array}
$$

    </div>
  </div>
</div>


\\
Where the left is for directly using a square root and right for a reciprocal sqrt.  So, anyway, just rewritting the conversion from the last post:

<div class="alert alert-warning" role="alert" markdown="1">
<b>WARNING:</b> this is unlikely to be what you want. It introduces extra computational complexity if the *swing* and *twist* aren't modified away from the axis and plane of the decomposition.
</div>

{% highlight c %}

// factor Q = Qr*Qz where Qz is a rotation about Z and Qr
// about an axis in the XY plane.
void quat_factor_z(quat_t* qz, quat_t* qr, quat_t* q)
{
  // cut-off for approaching Qz = identity & Qr
  // angle's magnitude is approaching pi.
  static const float THRESHOLD = (ULP1*ULP1);
  
  float x=q->x, y=q->y, z=q->z, w=q->w;
  float t=w*w+z*z;

  if (t > THRESHOLD) {
    float s  = 1.f/sqrtf(t);
    float fx = s*(w*x-y*z);
    float fy = s*(w*y+x*z);

    quat_set(qz,0.f,0.f, s*z, s*w);
    quat_set(qr,fx ,fy,  0.f, s*t);
  }
  else {
    quat_set_id(qz);
    quat_dup(qr,q);
  }
}

// reconstruct
void quat_recon_z(quat_t* q, quat_t* qz, quat_t* qr)
{
  quat_mul(q,qr,qz);
}

{% endhighlight %}

\\
And for what it's worth the probability of hitting the special case side of the branch is significantly less than (SEE: *"Cumulative volume of ..."* [here](https://marc-b-reynolds.github.io/quaternions/2017/11/10/AveRandomRot.html)):

$$ 1-\tfrac{1}{\pi} \left(2~\operatorname{acos}(T)-\operatorname{sin}(2~\operatorname{acos}(T)) \right)$$

where $T$ is the square root of `THRESHOLD` and it's significantly less because the equation is for when $w$ alone is approaching zero.  But anyway the code is using $T=2^{-24}$ and plugging that in gives $\approx 7.6 \times 10^{-8}$

<br>

------

Derivation
------

\\
Repeating the derivation from *Hopf* for *twist* then *swing* $\left(Q=Q_sQ_t\right)$ about some axis $\mathbf{a}$ (unit):

$$ 
\begin{align}
\mathbf{b} &= Q\mathbf{a}Q^{*}                         \label{e1} \\
Q_s        &= \left(\mathbf{ba}^*\right)^{\frac{1}{2}} \label{e2} \\
Q_t        &= \left(Q_s^*\right)Q                      \label{e3}
\end{align}
$$

* $\eqref{e1}$ finds how $Q$ rotates the axis $\mathbf{a}$.
* $\eqref{e2}$ finds the torque minimal (smallest magnitude angle) rotation that likewise maps $\mathbf{a}$ to $\mathbf{b}$. This is the *swing*. <br>This is solving: $b=Q_saQ_s^*$ for $Q_s$. (SEE: [here](https://marc-b-reynolds.github.io/quaternions/2016/08/09/TwoNormToRot.html) or from a vector perspective [here](https://marc-b-reynolds.github.io/quaternions/2018/06/12/SolveQRotVect.html))
* $\eqref{e3}$ finds the *twist*. Since $Q=Q_sQ_t$ multiply the left sides of both by $Q_s^*$ to solve for $Q_t$

and for the reversed ordered variant $\left(Q=Q_sQ_t\right)$ is found by $Q_r = Q\left(Q_s^*\right)$.

The twist can likewise be derived using plane reflections:

$$ 
\begin{align}
Q_s        &= \lVert \tfrac{1}{2}\left(Q-\mathbf{a}Q\mathbf{a}\right) \rVert \label{r2} \\
\end{align}
$$

* Reflect the bivector part by the $\mathbf{a}$. This is: $\mathbf{a}Q\mathbf{a}$
* isloate the part in $\mathbf{a}$ which is: $\frac{1}{2}\left(Q-\mathbf{a}Q\mathbf{a}\right)$
* normalize and done

and with the twist we can solve for the swing. I find this less satisfing but it's pretty easy to remember. Just zero out the bivector part independent of the desired axis and normalize.


<br>

------

Reconstruction
------

\\
To restore the quaternion we just need to perform the product $\eqref{def}$. 

Renaming each of the axis align factorizations:

$$ 
\begin{array}{lrl}
Q_{tx} = t_c + \bvec{t_s, 0, 0} & & Q_{sx} = s_c + \bvec{0, s_x, s_y} \\
Q_{ty} = t_c + \bvec{0, t_s, 0} & & Q_{sy} = s_c + \bvec{s_x, 0, s_y} \\
Q_{tz} = t_c + \bvec{0, 0, t_s} & & Q_{sz} = s_c + \bvec{s_x, s_y, 0}
\end{array}
$$

\\
Given these definations the product expand to:

$$ 
\begin{array}{lcl}
Q &= Q_{sx} Q_{tx} &= s_c t_c + \bvec{s_c t_s,           s_x t_c {\color{foo}+} s_y t_s, s_y t_c {\color{foo}-} s_x t_s} \\
  &= Q_{sy} Q_{ty} &= s_c t_c + \bvec{s_x t_c {\color{foo}-} s_y t_s, s_c t_s,           s_y t_c {\color{foo}+} s_x t_s} \\
  &= Q_{sz} Q_{tz} &= s_c t_c + \bvec{s_x t_c {\color{foo}+} s_y t_s, s_y t_c {\color{foo}-} s_x t_s, s_c t_s}
\end{array}
$$


\\
and the reversed order variants:

$$ 
\begin{array}{lcl}
Q &= Q_{tx} Q_{sx} &= s_c t_c + \bvec{s_c t_s,           s_x t_c {\color{foo}-} s_y t_s, s_y t_c {\color{foo}+} s_x t_s} \\
  &= Q_{ty} Q_{sy} &= s_c t_c + \bvec{s_x t_c {\color{foo}+} s_y t_s, s_c t_s,           s_y t_c {\color{foo}-} s_x t_s} \\
  &= Q_{tz} Q_{sz} &= s_c t_c + \bvec{s_x t_c {\color{foo}-} s_y t_s, s_y t_c {\color{foo}+} s_x t_s, s_c t_s}
\end{array}
$$

<br>

------

Comments on implementations
------

\\
Let's look at a messy strawman scalar version for filtering with everything explictly expanded. Something like this might be reasonable if there's only one decomp in the codebase but the point is the show the moving parts together.

{% highlight c %}

// helper funcs
inline uint32_t f32_to_bits(float x)
{
  uint32_t u; memcpy(&u, &x, 4); return u;
}

inline float f32_from_bits(uint32_t x)
{
  float f; memcpy(&f, &x, 4); return f;
}

// example decomposition with respect to X and
// reconstruction with some optional bell & whistles.

void swing_twist_filter_example(quat_t* q)
{
  static const bool normalize_w = true;
  static const bool restore_w   = true;

  static const float THRESHOLD = 0x1.0p-24f;

  //**********************************************
  // 'q' is the incoming quaternion to be filtered
  //**********************************************

  // just to make less visual noise below
  float x=q->x, y=q->y, z=q->z, w=q->w;

  // common sub-expression and test value (we really want math errno disabled!)
  float t=sqrtf(w*w+x*x);      // (*) change: x*x for other twist axis
  
  float ts,tc;                 // twist parameters: sin,cos about 'X' (unit complex number)
  float sx,sy,sc;              // swing parameters: sin*(px,py) & cos
                               //   where (px,py) is unit vector in plane ortho to 'X'
							   
  uint32_t sw = 0;             // isolated sign bit of 'w' (copysign removal)

  if (t > THRESHOLD) {
    float s0 = 1.f/t;

    // This is equivalent to negating Q if the scalar component (w) is negative
    // Only the twist parameters differ. The idea is for less complex filtering.
    if (normalize_w) {
      s1 = f32_from_bits(f32_to_bits(s0) | sw);
    //s1 = copysignf(s0,w);
    }
    
    // perform the decomposition
    sx = s0*(w*y-x*z);         // + for reversed order (*) axis dependent
    sy = s0*(w*z+x*y);         // - for reversed order (*) axis dependent
    sc =    t;
    ts = s1*x;
    tc = s1*w;
  }
  else {
    // limit case handling
    sx = y; sy = z; sc = 0;    // (*) axis dependent (sx,sy)
    ts = 0; tc = 1.f;
  }

  //**********************************************
  // FILTER GOES HERE
  // ...
  // END OF FILTERING
  //**********************************************

  // If 'w' has been normalized to positive it might be necessary to
  // restore the sign. An example is the quaternion could be a part of
  // a rotation sequence whos angle is greater than Pi to avoid an impossible
  // state.

    if (normalize_w && restore_w) {
      tc = f32_from_bits(f32_to_bits(tc) | sw);
      ts = f32_from_bits(f32_to_bits(ts) ^ sw);
    //tc  = copysignf(tc,w);
    //ts *= copysignf(1.f,w);
  }

  // recompose the filtered quaternion
  r->x = sc*ts;
  r->y = sx*tc + sy*ts;        // - for reversed order
  r->z = sy*tc - sx*ts;        // + for reversed order
  r->w = sc*tc;
}

{% endhighlight %}

\\
The commented out `copysignf` parts are due to current compliers (clang 13.0.0 & gcc 11.2) missing some reductions on these expressions so I just converted to manually performing the operations.

For codebases that need more than one factorization and/or libraries we'd want to do the obvious thing of having different functions for the decomposition and reconstuction. Notice that the various twist axis variants boil down to just variable renaming of each other so only one needs to be directly implemented and the others can just 'call' that version. Assuming a mildly agressive inlining the compiler will fully expand these.  A toy [header file](https://github.com/Marc-B-Reynolds/Stand-alone-junk/blob/master/src/SFH/swing_twist.h) does this.

Rewritting the above (including compile time configuation) gives:

{% highlight c %}
{
  tmtwist_t st;
  uint32_t  sign_w = 0;
  
  if (normalize_w)
    quat_to_swing_twist_x_n(&st, q, &sign_w);
  else
    quat_to_swing_twist_x(&st, q);
  
  // filtering here
  
  if (normalize_w && restore_w)
    swing_twist_x_to_quat_n(q, &st, sign_w);
  else
    swing_twist_x_to_quat(q, &st);
}

{% endhighlight %}
