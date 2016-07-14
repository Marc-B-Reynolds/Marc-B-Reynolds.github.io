---
layout:       post
title:        Orthonormal basis from normal via quaternion similarity
tagline:      
categories:   [quaternions]
tags:         [graphics, interpolation, quantization, distribution]
description:  two methods of computing an orthonormal basis from a unit (bi)vector.
---

<iframe width="600" height="400" src="http://www.shadertoy.com/embed/lldGRM" frameborder="0" allowfullscreen></iframe>{: .center-image }

\\
This is a quick note about some methods for computing an orthonormal basis from a unit (bi)vector.  The math directly follows that in "[Implied normals (unit bivectors as rotations)]({{site.base}}/quaternions/2016/06/26/QuatNormal.html)". 

\\
A re-cap of the needed bits from *implied*. Start with a standard orthonormal basis set:

$$
\begin{eqnarray*}
  \mathbf{x} & = & \left(1,0,0\right) \\ 
  \mathbf{y} & = & \left(0,1,0\right) \\ 
  \mathbf{z} & = & \left(0,0,1\right) 
\end{eqnarray*}
$$

\\
Given a unit (bi)vector:

$$\mathbf{v} = \left(v_x,~v_y,~v_z\right)$$

\\
Then the quaternion $Q$ which rotates $\mathbf{z}$ to $\mathbf{v}$:

$$
\begin{eqnarray}
  Q & = & \sqrt{\mathbf{v}\mathbf{z}^*}      \nonumber \\
    & = & \sqrt{-\mathbf{vz}}                \nonumber \\
    & = & \sqrt{\mathbf{z} \cdot \mathbf{v} + \mathbf{z} \times \mathbf{v}}      \nonumber \\
    & = & \sqrt{v_z + \left(-v_y,~v_x,~0 \right)}  \nonumber \\
    & = & \frac{1+v_z}{\sqrt{2+2v_z}} + \left(\frac{-v_y}{\sqrt{2+2v_z}},~\frac{v_x}{\sqrt{2+2v_z}},~0\right) \label{fx}
\end{eqnarray}
$$

<br>

------

Single reference direction <small>sticking with *up* seems reasonable</small>
------

\\
Directly using $\eqref{fx}$ we can transform (rotate) the standard orthonormal set {$\mathbf{x,y,z}$} into a new set {$\mathbf{x',y',v}$}:



$$
\begin{eqnarray}
  \mathbf{x}'  & = & Q\mathbf{x}Q^*  \nonumber \\
               & = & \left(  v_z+\frac{v_y^2}{v_z+1},~-\frac{v_x v_y}{v_z+1},~-v_x  \right) \label{xp} \\
                                     \nonumber \\
  \mathbf{y}'  & = & Q\mathbf{y}Q^*  \nonumber \\
               & = & \left(  -\frac{v_x v_y}{v_z+1},~1-\frac{v_y^2}{v_z+1},~-v_y  \right)   \label{yp} \\
                                     \nonumber \\
  \mathbf{z}'  & = & Q\mathbf{z}Q^* = \mathbf{v} \nonumber \\
	       
\end{eqnarray}
$$

\\
A simple translation of $\eqref{xp}$ and $\eqref{yp}$ into 'C' gives:

<br>
{% highlight c %}
void ortho_basis_1(vec3_t* v, vec3_t* xp, vec3_t* yp)
{
  float x = -v->x;
  float y =  v->y;
  float z =  v->z; 
  float a = y/(z+1.f);  // y/(z+1)
  float b = y*a;        // y^2/(z+1)
  float c = x*a;        // -xy/(z+1)
  
  vec3_set(xp, z+b, c,      x);  // {z+y/(z+1),   -xy/(z+1), -x}
  vec3_set(yp, c,   1.f-b, -y);  // {-xy/(z+1), 1-y^2/(z+1), -y}
}
{% endhighlight %}

\\
Simply inspecting the code or recalling the math shows a problem as $\mathbf{v}$ approaches $-\mathbf{z}$, we have an infinite number of solutions for $Q$ and numerically the math explodes. In cases were this is possible the above can be augmented by adding a branch and setting to `xp` and `yp` to two orthogonal vectors in the {$\mathbf{xy}$}-plane or the "correct" values in the limit.  Choice of branch point and values depend on use requirements. 

\\
The math basis for this method when using positive $\mathbf{z}$ as a single reference direction is identically to that in Frisvad[^frisvad].  However the math is presented differently and the derivation is carried through differently as well:

$$
\begin{eqnarray}
  \mathbf{x_f}  & = & \left(  1+\frac{v_x^2}{v_z+1},~-\frac{v_x v_y}{v_z+1},~-v_x  \right) \nonumber \\
  \mathbf{y_f}  & = & \left(  -\frac{v_x v_y}{v_z+1},~1-\frac{v_y^2}{v_z+1},~-v_y  \right) \nonumber 
\end{eqnarray}
$$

\\
which directly translation into code results at least four instead of two products.


<br>

------

Two opposite reference directions <small>*up* and *down*, why not?</small>
------

\\
A second possible solution is to form the rotation with nearest direction of two choices.  As an example we can form the rotation that maps $-\mathbf{z}$ to $\mathbf{v}$ and use that when $v_z < 0$:

$$
\begin{eqnarray}
  Q_n & = & \sqrt{\mathbf{vz}}              \nonumber \\
    & = & \frac{1-v_z}{\sqrt{2-2v_z}} + \left(\frac{v_y}{\sqrt{2-2v_z}},~-\frac{v_x}{\sqrt{2-2v_z}},~0\right) \label{nx}
\end{eqnarray}
$$

\\
Rotating $\mathbf{x}$ and $\mathbf{y}$ by $Q_n$:

$$
\begin{eqnarray}
  \mathbf{x}'  & = & \left( -v_z-\frac{v_y^2}{v_z-1},~\frac{v_x v_y}{v_z-1},~v_x  \right) \label{xn} \\
  \mathbf{y}'  & = & \left( \frac{v_x v_y}{v_z-1},~1+\frac{v_y^2}{v_z-1},~v_y  \right)   \label{yn}
\end{eqnarray}
$$

\\
The differences between $\eqref{xp} \eqref{yp}$ and $\eqref{xn} \eqref{yn}$ can be expressed by taking the absolute value of $v_z$ and products of $s_z= \text{sgn}\left(v_z\right)$. So we can combine the two expressions as follows: 

$$
\begin{eqnarray}
  \mathbf{x}'  & = & \left( \abs{v_z}+\frac{v_y^2}{\abs{v_z}+1},~-\frac{v_x v_y}{\abs{v_z}+1},~-s_z v_x \right) \label{xpn} \\
  \mathbf{y}'  & = & \left(-\frac{v_x v_y}{\abs{v_z}+1},~1-\frac{v_y^2}{\abs{v_z}+1},~-s_z v_y  \right) \label{ypn}  
\end{eqnarray}
$$

<br>

{% highlight c %}
static inline void ortho_basis_2(vec3_t* v, vec3_t* xp, vec3_t* yp)
{
  float x  = v->x;
  float y  = v->y;
  float z  = v->z;
  float sz = z >= 0.f ? -1.f : +1.f;  // -sgn(z) = -sz
  float az = fabsf(z);
  float a  = y/(az+1.f); //   y/(|z|+1)
  float b  = y*a;        // y^2/(|z|+1)
  float c  = -x*a;       // -xy/(|z|+1)
  
  vec3_set(xp, az+b, c,     sz*x); // {|z|+y/(|z|+1),   -xy/(|z|+1), -x}
  vec3_set(yp, c,    1.f-b, sz*y); // {  -xy/(|z|+1), 1-y^2/(|z|+1), -y}
}
{% endhighlight %}

\\
This method has a potential issue since $Q_n$ maps $-\mathbf{z}$ to $\mathbf{v}$.  This means the source basis set is {$\mathbf{x,y,-z}$} when $v_z < 0$ which is left-handed[^handed].  This is not an issue if we simply need two mutually orthogonal normals in the complement plane, however for cases where it is an issue we can simply negate one of the computed normals in the $v_z < 0$ case which is equivalent to multiplying through by $s_z$.  For both we then have:

$$
\begin{eqnarray}
  s_z \mathbf{x}'  & = & \left( v_z+\frac{s_z v_y^2}{\abs{v_z}+1},~-\frac{s_z v_x v_y}{\abs{v_z}+1},~-v_x \right) \label{xpnn} \\
  s_z \mathbf{y}'  & = & \left(-\frac{s_z v_x v_y}{\abs{v_z}+1},~s_z-\frac{s_z v_y^2}{\abs{v_z}+1},~-v_y  \right) \label{ypnn}  
\end{eqnarray}
$$

\\
Choosing to use $s_z \mathbf{x}'$ gives us:


{% highlight c %}
void ortho_basis_2a(vec3_t* v, vec3_t* xp, vec3_t* yp)
{
  float x  = -v->x;
  float y  = v->y;
  float z  = v->z; 
  float sz = z >= 0.f ? 1.f : -1.f;
  float a  = y/(fabsf(z)+1.f);
  float b  = y*a;
  float c  = x*a;
  
  vec3_set(xp, z+sz*b, sz*c, x); 
  vec3_set(yp, c,   1.f-b, -sz*y); 
}
{% endhighlight %}

<br>

------

Other methods
------

Thanks to Stephen Hill for providing references to some existing methods and (better yet) for writting up a summary in a blog post[^hill] and a link to a second pair of methods by Sam Hocevar[^lol].

These methods are all roughly:  Given $\mathbf{v}$ inspect the components and choose vector $\bar{\mathbf{w}}$ based on them.  Cross product creates an orthgonal vector.  This effectively partitions the sphere into "combing" groups.

The method above can be extended at increased computational complexity to shape how the normals are combed.  As an example we could add a twist in {$\\mathbf{x,y}$} like that found in a Wiechel projection or any of a number of craziness.

<br>

------

Viz
------

These figures roughly show the behavior of one output normal.  I'm not finding a quick-and-dirty viz that I really like.  Down positive z, positive x and negative z.

Method 1: Minimal angle orientation change on sphere.
![wtf]({{site.base}}/assets/figures/onbasis/method1.jpg 'WTF?'){: .center-image }

Method 2: Minimal angle orientation change on two half spheres.
![wtf]({{site.base}}/assets/figures/onbasis/method2.jpg 'WTF?'){: .center-image }

Hughes Moller: 
![wtf]({{site.base}}/assets/figures/onbasis/hughesmoller.jpg 'WTF?'){: .center-image }

------

References and Footnotes
------

[^handed]: If you don't like "handedness" then $\mathbf{x}' \times \mathbf{y}' = -\mathbf{v}$ instead of $\mathbf{v}$.
[^hm]:      *"Building an Orthonormal Basis from a Unit Vector"*, J. Hughes, T. Moller, Journal of Graphics Tools 4:4 (1999), 33-35.
[^stark]:  *"Efficient Construction of Perpendicular Vectors without Branching"*, M. Stark, Journal of Graphics Tools 14:1 (2009)
[^hill]:    *"Perpendicular Possibilities"*, Stephen Hill, 2011 [blog post](http://blog.selfshadow.com/2011/10/17/perp-vectors/)
[^frisvad]: *"Building an orthonormal basis from a 3d unit vector without normalization"*, Jeppe Frisvad, Journal of Graphics Tools 16:33, (2012) [PDF/code page](http://www.imm.dtu.dk/~jerf/papers/abstracts/onb.html)
[^lol]:      *"On picking an orthogonal vector (and combing coconuts)"*, Sam Hocevar, 2013 [blog post](http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts)