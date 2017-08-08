---
layout:       post
title:        Minimum magnitude angle rotation between two normals
tagline:      
categories:   [quaternions]
tags:         [rotations, normal, graphics]
description:  Given two normals find the rotation from one to the other
---

**Problem:** find the minimum magnitude angle rotation from $\hat{a}$ and $\hat{b}$, where both are unit vectors.

The unit quaternion $Q$ that rotates (unit bivectors) $a$ to $b$ can be expressed as:

$$
\mathbf{b} = Q~\mathbf{a}~Q^{-1}
$$

\\
and solving for $Q$ gives:

$$
\begin{eqnarray*}
  Q  &=& \sqrt{\mathbf{ba}^{-1}}  \\
     &=& \sqrt{\mathbf{ba}^*}  \\
     &=& \sqrt{-\mathbf{ba}}   \\
     &=& \sqrt{\mathbf{a} \cdot \mathbf{b} + \mathbf{a} \times \mathbf{b}} \\
     &=& \frac{1+\mathbf{a} \cdot \mathbf{b}}{\sqrt{2+2~\mathbf{a} \cdot \mathbf{b}}} + \frac{\mathbf{a} \times \mathbf{b}}{\sqrt{2+2~\mathbf{a} \cdot \mathbf{b}}} \\
     &=& \frac{1+d}{\sqrt{2+2d}} + \frac{\mathbf{v}}{\sqrt{2+2d}}
\end{eqnarray*}
$$

\\
This is well-covered material and if you've read my previous few posts then it is *copy-n-paste* math, *beat-a-dead-horse* material.  The bullet points are:

* The relative unit quaternion ($R$) between the unit bivectors is : $R=\mathbf{ba}^{-1}=...=\mathbf{a} \cdot \mathbf{b} + \mathbf{a} \times \mathbf{b}$. This maps the plane that spans $\mathbf{a},\mathbf{b}$ and 0 to the complex plane that spans ${1,R}$ such that $\mathbf{a}$ is mapped to 1.
* All that is left is taking the square root (half-angle) of $R$, which is a complex function so the remainder of the derivation is simply over complex numbers.
* When $b=-a$ we have $R=-1$ and $Q=\sqrt{-1}$. The math explodes since there are an infinite number of solutions instead of one unique (all unit bivectors).   Numerically the math explodes when approaching.

\\
Ignoring the degenerate case this translates into:

{% highlight glsl %}
vec4 q_from_normals(vec3 a, vec3 b)
{
  float k = 1.0+dot(a,b);         // 1+d
  float s = inversesqrt(k+k);     // 1/sqrt(2+2d)
  return vec4(s*cross(a,b), k*s); // (1+d)/sqrt(2+2d) + (a x b)/sqrt(2+2d)
}

{% endhighlight %}

<br>

Performing some renaming of variables:

$$
\mathbf{v}=\left(x,y,z\right) \\
r=\frac{1}{1+d}
$$

and applying $Q$ to the standard orthonormal basis yields a rotation matrix:

$$
{
\left( \begin{array}{ccc}
1-r(y^2+z^2) & rxy+z        & rxz-y \\
rxy-z        & 1-r(x^2+z^2) & ryz+x \\
rxz+y        & ryz-x        & 1-r(x^2+y^2)
\end{array} \right)
}
$$

\\
The above is reduced in what I consider the most obvious way, but the diagonals have that horrible one minus.  Replacing the one with $Q \cdot Q$ and reducing yields the much nicer:

$$
{
\left( \begin{array}{ccc}
d+rx^2 & rxy+z  & rxz-y \\
rxy-z  & d+ry^2 & ryz+x \\
rxz+y  & ryz-x  & d+rz^2
\end{array} \right)
}
$$

\\
Since we can rewrite $r$ as:

$$ \frac{1}{1+d} = \frac{1-d}{1-d^2}$$

\\
shows that the above is identical to the non-degenerate case presented by Moller and Hughes[^mh99].

Some alternate constructions for the same material:

* Inigo Quilez walks to the same matrix solution as an example of the dangers of composing functions instead of carrying through derivations.[^iq]
* Sam Hocevar builds the quaternion solution from vector calculus and trigonometry.[^sh]

Toy C code can be found: [HERE](http://github.com/Marc-B-Reynolds/Stand-alone-junk/blob/master/src/Posts/normals_to_rot.c)

------

References and Footnotes
------

[^mh99]:   *"Efficiently Building a Matrix to Rotate One Vector to Another"*, Tomas Moller, John F. Hughes, 1999, ([PDF](http://cs.brown.edu/~jfh/papers/Moller-EBA-1999/paper.pdf))
[^iq]:     *"avoiding trigonometry"*, Inigo Quilez, 2013, ([PAGE](http://www.iquilezles.org/www/articles/noacos/noacos.htm))
[^sh]: *"Beautiful maths simplification: quaternion from two vectors"*, Sam Hocevar, 2013 ([PAGE](http://lolengine.net/blog/2013/09/18/beautiful-maths-quaternion-from-vectors))
