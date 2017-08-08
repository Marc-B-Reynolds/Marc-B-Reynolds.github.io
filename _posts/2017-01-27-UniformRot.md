---
layout:       post
title:        Uniform random 3D rotations
tagline:      
categories:   [distribution]
tags:         [random, rotations]
description:  a brief note on generating random rotations in 3D
---

This is a quick extension to a previous post on uniform points on disc, circle, sphere and spherical caps[^uniform].  We can also move away from trig or the rejection method by approximing equal-area disc[^sdisc].  I'll assume basic understanding of representing rotations by quaternions.  For a matrix simply follow through the deriviation by applying the similarity transform to an orthogonal basis.

Informally there are three general types of random rotations in 3D...one, two and three degrees of freedom.

------

About fixed axis <small>one degree of freedom</small>
------
{:#1dof}

\\
The simplest general case is generating a random rotation about a fixed axis (considered as unit bivector $\mathbf{u}$) which is equivalent to the two dimensional case.  We can generate a point on the half-circle (positive in $x$) and set the scalar part of the quaternion to the $x$ coordinate and set the bivector part to $\mathbf{u}$ scaled by the $y$ result.
<br>

{% highlight c %}
void uniform_quat_about_u(quat_t* q, vec3_t* u)
{
  vec2_t p;
  float  d = uniform_circle_px(&p);
  quat_set(q, p.y*u->x, p.y*u->y, p.y*u->z, p.x);
}
{% endhighlight %}

\\
For non-uniform, if the distribution on the half-circle has angle probability function $P(\theta)$, then the result probability of the rotation angle is $P(2\theta)$.

<br>

------

With a fixed reference direction <small>two degrees of freedom</small>
------
{:#2dof}

\\
In axis-angle speak moving to two degress of freedom fixes the axis-of-rotation to a plane. We can choose a fixed reference direction which is orthogonal to that plane  (say positive $\mathbf{z}$), generate a uniform point on the sphere and find the rotation. 

First generate a uniform point on the disc: $p=\left(x,~y\right)$ and map to the sphere[^uniform]:

$$ \left( x,~y \right) \rightarrow \left(2x\sqrt{1-\left(x^2+y^2\right)},~2y\sqrt{1-\left(x^2+y^2\right)},~1-2\left(x^2+y^2\right)\right) $$

\\
Find the quaternion that rotates $\left(0,0,1\right)$ to some $\left(x,~y,~z\right)$[^inormals]:

$$ \left( x,~y,~z \right) \rightarrow \frac{1}{2+2z}\left(1+z + \left(-y,~x,~0\right)\right) $$

\\
Composing the previous operations results in:

$$ \left( x,~y \right) \rightarrow \sqrt{1-\left(x^2+y^2\right)} + (-y, ~x, 0) $$

Breaking this down we have:

* the direction of the axis of rotation:  $(-y,~x,0)$
  * CCW orthogonal to point on disc (which is logically in $\mathbf{xy}$-plane)
  * magnitude is logically $\sin\left(\frac{\theta}{2}\right)$
* the scalar part is logically $\cos\left(\frac{\theta}{2}\right) = \sqrt{1-\sin(\theta)^2} = \sqrt{1-d}$

We could have directly generated this result but requires figuring out how to set up a Jacobian...and I hate doing that. Since our values on the disc are random, we can also simply reorder the bivector part and drop the sign (not rotate) without any impact, giving us:

$$ \left( x,~y \right) \rightarrow \sqrt{1-\left(x^2+y^2\right)} + (x,y,0) $$

{% highlight c %}
void uniform_quat_from_z(quat_t* q)
{
  vec2_t p;
  float  d = uniform_disc(&p);
  float  s = sqrtf(1-d);            // root on (0,1] can approx + fix-up
  quat_set(q, p.x, p.y, 0.f s);
}
{% endhighlight %}

\\
For non-uniform (with the reordering and sign drop) then if the probability on the disc is $P(x,y)$ is the probability of a rotation of $2 \text{ acos}\left(1-(x^2+y^2)\right)$ about axis in direction $(x,y,0)$. More generally we can shape the disc distribution and convert to cap.

<br>

------

The Real Deal (arbitrary) <small>three degrees of freedom</small>
------
{:#3dof}

\\
George Marsaglia in same paper[^gm] as uniform points on the 3D sphere $\left(\mathbb{S}^2\right)$ hand-wavingly[^pervognsen] presented a method for the 4D sphere $\left(\mathbb{S}^3\right)$.  

Given two uniform points on the unit disc: $~p_0,~p_1$ and their norm: $d_n=p_n \cdot p_n$:


$$
   \left( x_0,~y_0,
   ~x_1 \sqrt{\frac{1-d_0}{d_1}},
   ~y_1 \sqrt{\frac{1-d_0}{d_1}}
   \right)
$$

\\
Making an arbitrary choice for the scalar part ($x_0$) the quaternion is:

$$
   x_0 + \left(y_0,
   ~x_1 \sqrt{\frac{1-d_0}{d_1}},
   ~y_1 \sqrt{\frac{1-d_0}{d_1}}
   \right)
$$

{% highlight c %}
void uniform_quat(quat_t* q)
{
  vec2_t p0,p1;
  float  d1 = uniform_disc(&p1) + EPS;
  float  s1 = rsqrt(d1);
  float  d0 = uniform_disc(&p0);  // or positive in 'x' since -Q & Q are equivalent
  float  s0 = sqrtf(1.f-d0);
  float  s  = s0*s1;

  quat_set(q, p0.y, s*p1.x, s*p1.y, p0.x);
}
{% endhighlight %}

<br>

------

References and Footnotes
------

[^gm]:       *"Choosing a Point from the Surface of a Sphere"*, George Marsaglia, 1972 ([PDF](http://projecteuclid.org/euclid.aoms/1177692644))
[^sdisc]:    *"Square/Disc mappings"*, ([local post]({{site.base}}/math/2017/01/08/SquareDisc.html))
[^inormals]: *"Implied normals (unit bivectors as rotations)"*, ([local post]({{site.base}}/quaternions/2016/06/26/QuatNormal.html))
[^uniform]:  *"Uniform points on disc, circle, sphere and caps"*, ([local post]({{site.base}}/distribution/2016/11/28/Uniform.html))
[^pervognsen]: Per Vognsen provides a simple justification in this gist: ([link](http://gist.github.com/pervognsen/32ec4841ef53346f65c5033c3bd262b6))
