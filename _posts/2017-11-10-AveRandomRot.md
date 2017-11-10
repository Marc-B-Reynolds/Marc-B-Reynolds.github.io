---
layout:       post
title:        Volume element of SO(3) and average uniform random rotation angle
categories:   [quaternions]
tags:         [rotation, quantization, distribution, sampling, map]
description:  Find the expression of a volume element of rotations.
---

\\
While toying around with some bits-n-bobs for my next quaternion quantization post I made a joke twitter post:

> Pointless 3D fact of the day:
> Given uniform random rotations, then their average is
> $\pi/2$ + $2/\pi$ radians (~126.476 degrees)

\\
I was probably type parts of this up anyway, so I'll just do it now...sorry it's kinda brain-dump like.

------

Prelim
------

\\
Rehashing the same old stuff.  The set of all quaternions $\mathbb{H}$, excluding zero, can represent 3D rotations by the similarity transform.  The subset of quaternions representing the same rotation form a line through the origin (excluded).  Limiting ourselves to unit quaternions this line intersects the 4D sphere (3-sphere or $\mathbb{S}^3$) at two points $Q$ and $-Q$.  By further limiting ourselves to unit quaternions with positive (or zero) scalars then each point in this subset represents a unique 3D rotation.  So unit quaternions lie on $\mathbb{S}^3$ and we can represent all rotations by half of this sphere.  By a web-search or hitting up a math reference site we can finding the surface volume (or 3-dimensional cubic hyperarea if that's your thing) of the sphere which is $2\pi^2$, so the volume of our half sphere is $\pi^2$.  Gotta run through the math of this since we need the volume element.  

EDIT:
* I should have checked [MathWorld](http://mathworld.wolfram.com/Hypersphere.html) for the volume element.
* And searched [math.stackexchange](https://math.stackexchange.com/questions/464419/mean-value-of-the-rotation-angle-is-126-5%C2%B0) for average rotation.

<br>

------

Volume of $SO\left(3\right)$
------

\\
Define a unit quaternion with the bivector part in azimuth $\left(\alpha\right)$ / inclination $\left(\beta\right)$ spherical coordinates:

$$ Q = \cos\left(\theta\right) + \sin\left(\theta\right) 
\Big(
 \sin\left(\alpha\right)\cos\left(\beta\right),
~\sin\left(\alpha\right)\sin\left(\beta\right), 
~\cos\left(\alpha\right) 
\Big) $$

\\
with angle ranges of:

$$ 
\theta  \in \left[0,~\frac{\pi}{2}\right],
~\alpha \in \left[0,~\pi\right],
~\beta  \in \left[0,~2\pi\right]
$$

\\
Now we have this 3D spherical coordinate thing that's be converted into a 4D Euclidean space thing and I want to compute the Jacobian determinant.  To get a square Jacobian (so we can compute the determinate) I'll just extend the quaternion to 4D by allowing non-unit magnitude (or a radius if you want) and multply the components by $r$:

$$ \begin{eqnarray*}
q_w & = & r~ \cos\left(\theta\right) \\
q_x & = & r~ \sin\left(\theta\right) \sin\left(\alpha\right)\cos\left(\beta\right) \\
q_y & = & r~ \sin\left(\theta\right) \sin\left(\alpha\right)\sin\left(\beta\right) \\
q_z & = & r~ \sin\left(\theta\right) \cos\left(\alpha\right)
\end{eqnarray*} $$

\\
Then we have the Jacobian:

$$
J = \left(
\begin{array}{cccc}
 \cos (\theta) & -r \sin (\theta ) & 0 & 0 \\
 \cos (\beta) \sin (\alpha ) \sin (\theta ) & r \cos (\beta ) \cos (\theta )
   \sin (\alpha) & r \cos (\alpha ) \cos (\beta ) \sin (\theta ) & -r \sin
   (\alpha ) \sin (\beta ) \sin (\theta ) \\
 \sin (\alpha ) \sin (\beta ) \sin (\theta ) & r \cos (\theta ) \sin (\alpha )
   \sin (\beta ) & r \cos (\alpha ) \sin (\beta ) \sin (\theta ) & r \cos (\beta
   ) \sin (\alpha ) \sin (\theta ) \\
 \cos (\alpha ) \sin (\theta ) & r \cos (\alpha ) \cos (\theta ) & -r \sin
   (\alpha ) \sin (\theta ) & 0 \\
\end{array}
\right)
$$

\\
Setting $r=1$ then the volume element is:

$$
dA = \left\vert \det J \right\vert d\theta~d\alpha~d\beta = \sin\left(\theta\right)^2~\sin\left(\alpha\right)~d\theta~d\alpha~d\beta
$$

\\
The volume is then:

$$
A = \int_{0}^{2\pi} \int_{0}^{\pi} \int_{0}^{\pi/2} \sin\left(\theta\right)^2~\sin\left(\alpha\right)~d\theta~d\alpha~d\beta = \pi^2
$$

\\
Great it agrees with what we could have looked up on Wikipedia.

<br>

------

Volume of $SO\left(3\right)$ <small>take 2</small>
------

\\
I don't need it but since I'm typing we could have done this a second way.  Without adding an extra dimension the Jacobian is:

$$
J = \left(
\begin{array}{ccc}
 -\sin (\theta ) & 0 & 0 \\
 \cos (\beta ) \cos (\theta ) \sin (\alpha ) & \cos (\alpha ) \cos (\beta ) \sin (\theta ) & -\sin
   (\alpha ) \sin (\beta ) \sin (\theta ) \\
 \cos (\theta ) \sin (\alpha ) \sin (\beta ) & \cos (\alpha ) \sin (\beta ) \sin (\theta ) & \cos
   (\beta ) \sin (\alpha ) \sin (\theta ) \\
 \cos (\alpha ) \cos (\theta ) & -\sin (\alpha ) \sin (\theta ) & 0 \\
\end{array}
\right)
$$

\\
then the Euclidean metric:

$$
g = J^T J = \left(
\begin{array}{ccc}
 1 & 0 & 0 \\
 0 & \sin ^2(\theta ) & 0 \\
 0 & 0 & \sin ^2(\alpha ) \sin ^2(\theta ) \\
\end{array}
\right)
$$

\\
and it's determinate:

$$ \det\left(g\right) = \sin\left(\theta\right)^4~\sin\left(\alpha\right)^2 $$

\\
The volume element is then:

$$
dA = \sqrt{\det\left(g\right)}~d\theta~d\alpha~d\beta = \sin\left(\theta\right)^2~\sin\left(\alpha\right)~d\theta~d\alpha~d\beta
$$

Same as before.

<br>

------

Average rotation angle of uniform random rotations
------

\\
Not explictly mentioned so far is that $\theta$ is the measure in quaternion space which translates into a $2\theta$ rotation. I want to find the average rotation angle $\phi$ of a uniform distribution of random rotations.  We have $q_w = \cos\left(\frac{\phi}{2}\right) = \cos\left(\theta\right)$ and the expected value is:

$$ \begin{eqnarray*}
E\lbrack\phi\rbrack & = & 2~E\lbrack\theta\rbrack     \\
                    & = & \frac{1}{A}\int 2~\theta~dA \\
                    & = & \frac{2}{\pi^2}\int ~\theta~\sin\left(\theta\right)^2~\sin\left(\alpha\right)~d\theta~d\alpha~d\beta \\
                    & = & \frac{4}{\pi^2} \left(\pi + \frac{\pi^3}{4}\right) \\
                    & = & \frac{\pi}{2} + \frac{2}{\pi}
\end{eqnarray*} $$

There's surely an easier way to do this and I'm shaky on the math here, so an empirical validation: [CLICK](https://gist.github.com/Marc-B-Reynolds/68ad708c950f57f0e38a445f9e9ef697)
