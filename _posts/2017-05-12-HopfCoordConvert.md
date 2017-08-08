---
layout: post
title:      'Hopf coordinate conversion and torque minimal factorization'
categories: [quaternions]
tags :      [hopf,rotations]
description : Gives a mile-high overview of the Hopf map from geometry perspective.
plotly: true
---

In a number of my posts I've been using properties related to Hopf coordinates and maps without explictly point it out.  This is partly because I'm lazy and partly because I don't know the math very well. Since I'm going to use Hopf coordinates soon I'll give an informal and inexact overview using my kinda language and notation.

At the same time it makes sense to show some practical implementations of conversion between quaternions and Hopf coordinates and factorizing a quaternion to a rotation pair: torque minimal and a twist.

For this post I'll be using angle measurements in $\mathbb{R}^3$ instead of my normally preferred in $\mathbb{H}$.

------

Hopf map
------

\\
We can represent rotations by unit quaternions which we can associate with the 4D unit sphere $ \left( \mathbb{S}^3 \right) $.  The Hopf map $\left(h\right)$ sends a point on $\mathbb{S}^3 $ to a point on the 3D unit sphere $\left(\mathbb{S}^2\right)$.

This section broken down into two parts.  The first is just a bullet-point list using geometry and rotations.  The second covers the same ground using algebra for brushstrokes of justification. 
$$ \newcommand{\set}[1]{\left\{ #1 \right\}} $$

### Rotation brushstrokes

* We choose a reference direction and the point on the unit sphere in that direction. I'm using $\mathbf{z}=\left(0,0,1\right)$.
* All directions orthogonal to the reference are in its orthogonal complement.  In this case the $\set{\mathbf{xy}}$-plane.
* Given a 3D rotation represented as unit quaternion $Q$ the rotation can be factored into a pair of rotations. One about the axis in our reference direction $\left(Q_t\right)$, the other about an axis in the orthogonal complement $\left(Q_r\right)$ which compose as $Q=Q_rQ_t$.
* The Hopf map is simply the point returned by rotating our reference by $Q$ and the result is independent of the $Q_t$ factor. 
  * The map sends an infinite number of rotations (all the possible $Q_t$) to each point on the 3D sphere ($p$) and the collection of all them is called a *fiber*.
  * All rotations about $\mathbf{z}$ are mapped to $\mathbf{z}$.  The torque minimal $Q_r$ factor is one.
  * All maximum magnitude rotations are mapped to $-\mathbf{z}$.  The torque minimal $Q_t$ factor is one.
  * Generally the value of $h\left(Q\right)=\mathbf{p}$ is uniquely determined by the $Q_r$ factor.

Applying the factors as a global rotation sequence:  Rotating $\mathbf{z}$ by any rotation about $\mathbf{z}$ leaves it in place.  Then rotating about any axis in the $\set{\mathbf{xy}}$-plane sends it to $\mathbf{p}$.


<br>

### Algebra brushstrokes

\\
The Hopf map $h: \mathbb{S}^3 \rightarrow \mathbb{S}^2$ can be directly formulated by choosing a reference direction and apply the similarity transform to find the point $\mathbf{p}$.  Or more simply: we treat $Q$ as a rotation, use it to rotate $\mathbf{z}$ and the result of the rotation is $\mathbf{p}$:

\\
Given a unit quaternion:

$$ \begin{equation} \label{Qr}
Q = w + \left(x,~y,~z\right) 
\end{equation} $$ 

\\
then the Hopf map can be expressed as:

$$
\begin{align}
h\left(Q\right) & =  Q\mathbf{z}Q^{*}  \nonumber \\
           & =  \left( 2(wy+xz), ~2(yz-wx), ~(w^2+z^2)-(x^2+y^2) \right) \nonumber \\
           & =  \left( 2(wy+xz), ~2(yz-wx), ~1-2(x^2+y^2) \right) \label{Pr} \\
           & =  \mathbf{p} \nonumber
\end{align}

$$

\\
The *fiber* of $\mathbf{p}$ is the set of all unit quaternions that map to it.  Let's call a rotation about the axis in the reference direction a *twist* and given the choice of $\mathbf{z}$ we can represent a *twist* of $\gamma$ as:

$$ \begin{equation} \label{rotZ}
Q_t =  \cos\left(\frac{\gamma}{2}\right)+\sin\left(\frac{\gamma}{2}\right)\mathbf{z}
\end{equation} $$ 

\\
where the angle $\gamma$ uniformally parameterizes the 2D unit circle $\left(\mathbb{S}^1\right)$ in the plane spanned by $ \set{ 1,\mathbf{z} } $.  Plugging into the map yields:

$$
\begin{align*}
h\left(Q_t\right)
 & = \left(0,~0, \cos\left(\frac{\gamma}{2}\right)^2+\sin\left(\frac{\gamma}{2}\right)^2 \right) \\
 & = \left(0,~0,~1\right) \\
 & = \mathbf{z}
\end{align*}
$$

\\
A straight forward result from our formulation of $h$.

Given some $\mathbf{p}=\left(a,~b,~c\right)$ we can find the torque minimal rotation $\left(Q_r\right)$ that maps (rotates) $\mathbf{z}$ to $\mathbf{p}$:

$$ 
\begin{align}
Q_r &= \sqrt{\mathbf{pz}^*} \nonumber \\
    &= \frac{1+c}{\sqrt{2(1+c)}} + \left(\frac{-b}{\sqrt{2(1+c)}}, ~\frac{a}{\sqrt{2(1+c)}}, ~0 \right) \label{rotM}
\end{align}
$$

\\
With the $p$ from $\eqref{Pr}$:

$$ \begin{equation} \label{rotQm}
Q_r = \sqrt{w^2+z^2} + \left(\frac{wx-yz}{\sqrt{w^2+z^2}},~\frac{wy+xz}{\sqrt{w^2+z^2}},~0\right) \\
\end{equation} $$

\\
And we can rewrite in terms of a rotation of $\alpha$ about an axis in the $\set{\mathbf{xy}}$-plane parameterized by $\beta$ (awkward choice here is to logically match $\eqref{rotM}$):

$$ \begin{equation} \label{rotZp}
Q_r =  \cos\left(\frac{\alpha}{2}\right)+\sin\left(\frac{\alpha}{2}\right)\left(-\sin\left(\beta\right),~\cos\left(\beta\right),~0\right)
\end{equation} $$

\\
So we have a unit quaternion which represents a rotation about an axis in the $\set{\mathbf{xy}}$-plane which we can associate with a point on the 3D unit sphere $\left(\mathbb{S}^2\right)$.

\\
The fiber of $\mathbf{p}$ can then be expressed as the composition of the two rotations $\eqref{rotZ}$ expanded with $\eqref{rotM}$ and $\eqref{rotZp}$:

$$
\begin{align*}
Q &= Q_r Q_t \\
  &= \frac{1}{2\left(1+c\right)} \left(\left(1+c\right) \cos\left(\frac{\gamma}{2}\right) + 
  \left(
    a \sin\left(\frac{\gamma}{2}\right) - b \cos\left(\frac{\gamma}{2}\right),
    a \cos\left(\frac{\gamma}{2}\right) + b \sin\left(\frac{\gamma}{2}\right),
  \left(1+c\right) \sin\left(\frac{\gamma}{2}\right)\right)\right) \\
  \\
  &= \cos\left(\frac{\alpha}{2}\right)\cos\left(\frac{\gamma}{2}\right) + \left(
  \sin\left(\frac{\alpha}{2}\right)\sin\left(\frac{\gamma}{2}-\beta\right),
  \sin\left(\frac{\alpha}{2}\right)\cos\left(\frac{\gamma}{2}-\beta\right),
  \cos\left(\frac{\alpha}{2}\right)\sin\left(\frac{\gamma}{2}\right)
\right)

\end{align*}
$$

\\
Simply plugging in the composed rotation into the map gives:

$$
\begin{align*}
h\left(Q\right) 
&=  Q\mathbf{z}Q^{*}  \\
&=  Q_r Q_t \mathbf{z} \left(Q_r Q_t\right)^{*}  \\
&=  Q_r Q_t \mathbf{z} Q_t^* Q_r^* \\
&=  Q_r \mathbf{z} Q_r^* \\
&=  \mathbf{p}
\end{align*}
$$

\\
So the result of $h\left(Q\right)$ is determined by the $Q_r$ factor.

<br>

------

Hopf coordinates
------

\\
Given the above we can represent a rotation a pair of coordinates, one on the unit sphere and another on a unit circle.  These coordinates parameterize two rotations which when composed represent the set of all 3D rotations can be parametrized by a triple of angles $\left(\alpha,\beta,\gamma\right)$:


$$
\begin{align*}
\gamma &= 2 \text{ atan} \left(\frac{z}{w} \right) \\
\alpha &= 2 \text{ asin} \left(\sqrt{x^2+y^2} \right) \\
       &= 2 \text{ acos} \left(\sqrt{w^2+z^2} \right) \\
       &= 2 \text{ atan} \left(\sqrt{\frac{x^2+y^2}{w^2+z^2}} \right) \\
\beta  &=   \frac{\gamma}{2} - \text{ atan} \left(\frac{x}{y}\right) \\
       &=   \text{ atan} \left(\frac{yz-wx}{wy+xz}\right)
\end{align*}
$$

A performance sub-optimal scalar version:

{% highlight c %}
void quat_to_hopf(vec3_t* v, quat_t* q)
{
  float x=q->x, y=q->y, z=q->z, w=q->w;

  // both gamma and alpha can be reduced to single parameter
  // atan by manually performing some argument reduction.
  // Also using a select to logically compute alpha via 
  // asin or acos drops a divide (either way one sqrt is
  // sufficient) by choosing the larger input parameter.
  // Okay beta can be as well...point being the other two
  // are already have less cases to be handled.
  float gamma = 2.f*atan2f(z, w);
  float alpha = 2.f*atan2f(sqrtf(x*x+y*y), sqrtf(w*w+z*z));
  float beta  = atan2f(y*z-w*x, w*y+x*z);

  v->x = alpha;
  v->y = beta;
  v->z = gamma;
}

void hopf_to_quat(quat_t* q, vec3_t* v)
{
  float ha    = 0.5f*v->x;
  float hg    = 0.5f*v->z;
  float b     = v->y;
  float sinHa = sinf(ha);
  float cosHa = cosf(ha);
  
  q->x = sinHa*sinf(hg-b);
  q->y = sinHa*cosf(hg-b);
  q->z = cosHa*sinf(hg);
  q->w = cosHa*cosf(hg);
}

{% endhighlight %}




<br>

------

Torque minimal factorization
------

\\
Given that an input rotation $Q$ can be factor into:

$$
Q = Q_r Q_t
$$

\\
We can solve for $Q_t$ to complete the factorization:

$$
\begin{align*}
Q_t &= \left(Q_{m}\right)^*Q \\
    &= \frac{w}{\sqrt{w^2+z^2}} + \left(0,~0,~\frac{z}{\sqrt{w^2+z^2}}\right) 
\end{align*}
$$

\\
A possible application of such a factorization is to filter and/or cap the *twist* for a sequence of rotations.  Toy scalar code example of the factorization:


{% highlight c %}

typedef struct {
  float mx,my,ms;   // wrt reference direction
  float tz,ts;      // twist
} tmtwist_t;

#define THRESHOLD (ULP1*ULP1)

void quat_to_tmtwist(tmtwist_t* d, quat_t* q)
{
  float x=q->x, y=q->y, z=q->z, w=q->w;
  float t=w*w+z*z;

  // For the full sphere the numerically unstable
  // portion could be improved.
  if (t > THRESHOLD) {
    float s = 1.f/sqrtf(t);
    d->mx = s*(w*x-y*z);
    d->my = s*(w*y+x*z);
    d->ms = s*t;
    d->tz = s*z;
    d->ts = s*w;
  }
  else {
    d->mx = x;
    d->my = y;
    d->ms = 0;
    d->tz = 0;
    d->ts = 1.f;
  }
}

void tmtwist_to_quat(quat_t* q, tmtwist_t* s)
{
  q->x = s->mx*s->ts + s->my*s->tz;
  q->y = s->my*s->ts - s->mx*s->tz;
  q->z = s->ms*s->tz;
  q->w = s->ms*s->ts;
}

{% endhighlight %}


<br>

------

References and Footnotes
------

[^pos]:     Otherwise replace $w$ by $\abs{w}$ and add $\text{sgn}\left(w\right)$ to the final scale value.  In this case the inverse transform returns $-Q$ when $w$ is negative.

<script>

// x-position for error traces
var xpos = new Array(100);
for (var i=0; i<100; i++) { xpos[i] = i*(1.0/99.0); };



</script>
