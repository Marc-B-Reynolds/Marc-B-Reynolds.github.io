---
layout: post
title:      "Factor a quaterion into a rotation about fixed axis and its orthogonal complement"
categories: [quaternions]
tags :      [rotations]
description : 'Given a unit quaternion $Q$ factor into a pair such that $Q = Q_rQ_a$ where $a$ is a predefined axis and $r$ is an orthogonal axis'
plotly: true
---

This is expanding on a detail from a previous post ([*"Hopf coordinate conversion and torque minimal factorization"*](https://marc-b-reynolds.github.io/quaternions/2017/05/12/HopfCoordConvert.html)) which is that a quaternion $Q$ can be factored into:

$$Q = Q_{xy}Q_z$$

\\
where $Q_z$ is a rotation about the $\mathbf{z}$-axis and $Q_{xy}$ is some orthogonal axis in the $\mathbf{XY}$-plane. There is example code for doing this (`quat_to_tmtwist`) but it is storing the results in a specialized structure that contains both factors and skips on the zero components. The inverse function (`tmtwist_to_quat`) is just a specialized product with this knowledge.  A couple of days ago I pointed someone to this info and on re-reading it was like...dude...not obvious.

So let's skim: given a unit quaternion:

$$ \begin{equation} \label{Qr}
Q = w + \left(x,~y,~z\right) 
\end{equation} $$ 

\\
we can factor into pair with respect to any rotation axis. For the three "coordinate" axes we have:

$$ 
\begin{align}
Q_x &= \frac{w}{\sqrt{w^2+x^2}} + \left(\frac{x}{\sqrt{w^2+x^2}},~0,~0\right) 
& Q_{yz} &= \sqrt{w^2+x^2} + \left(0,~\frac{wy-xz}{\sqrt{w^2+x^2}},~\frac{wz+xy}{\sqrt{w^2+x^2}}\right) \label{fx} \\

Q_y &= \frac{w}{\sqrt{w^2+y^2}} + \left(0,~\frac{y}{\sqrt{w^2+y^2},~0}\right) 
& Q_{xz} &= \sqrt{w^2+y^2} + \left(\frac{wx+yz}{\sqrt{w^2+y^2}},~0,~\frac{wz-xy}{\sqrt{w^2+y^2}}\right) \label{fy} \\

Q_z &= \frac{w}{\sqrt{w^2+z^2}} + \left(0,~0,~\frac{z}{\sqrt{w^2+z^2}}\right) 
& Q_{xy} &= \sqrt{w^2+z^2} + \left(\frac{wx-yz}{\sqrt{w^2+z^2}},~\frac{wy+xz}{\sqrt{w^2+z^2}},~0\right) \label{fz} \\

\end{align}
$$

\\
I can't see any reason to fully repeat the algebra (SEE: *Hopf*). The factorizations with respect to $\mathbf{x}\,\eqref{fx}$ and $\mathbf{y}\,\eqref{fy}$ simply are rinse and repeat. I'll stick with using $\mathbf{z}\,\eqref{fz}$ as the working example.

The factorization has the following properties: both $Q$ and $Q_{xy}$ rotate $\mathbf{z}$ (and the z component of a vector) to the same value:

$$ \mathbf{z}' = Q~\mathbf{z}~Q^* = Q_{xy}~\mathbf{z}~Q_{xy}^*$$

and $Q_z$ leaves $\mathbf{z}$ (and z components) unmodified:

$$  Q_z~\mathbf{z}~Q_z^* = \mathbf{z} $$

Rewritting the factors using two temporary variables we have:

$$ 
\begin{array}{rl}
t &=  w^2+z^2 \\
s &=  1/\sqrt{t} \\
Q_z &= sw + \left(0,~0,~sz\right) \\
Q_{xy} &= st + \left(\left(wx-yz\right),~\left(wy+xz\right),~0\right) \\

\end{array}
$$

\\
Now all of these factorizations have an ill conditioned region. When $t$ is approaching zero the input is approaching the maximum magnitude rotation ($\pi$) and axis is likewise approaching being independent of one one we're attempting to factor (again $\mathbf{z}$ for the working example).  So, anyway, just rewritting the conversion from the last post:

<br>

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

// just for show
void quat_recon_z(quat_t* q, quat_t* qz, quat_t* qr)
{
  quat_mul(q,qr,qz);
}

{% endhighlight %}

\\
And for what it's worth the probability of hitting the special case side of the branch is significantly less than (SEE: *"Cumulative volume of ..."* [here](https://marc-b-reynolds.github.io/quaternions/2017/11/10/AveRandomRot.html)):

$$ 1-\frac{1}{\pi} \left(2~\operatorname{acos}(T)-\operatorname{sin}(2~\operatorname{acos}(T)) \right)$$

where $T$ is the square root of `THRESHOLD` and it's significantly less because the equation is for when $w$ alone is approaching zero.  But anyway the code is using $T=2^{-24}$ and plugging that in gives $\approx 7.6 \times 10^{-8}$
