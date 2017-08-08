---
layout:       post
title:        Orthonormal basis from normal via quaternion similarity
categories:   [quaternions]
tags:         [graphics, interpolation, quantization, distribution]
description:  two methods of computing an orthonormal basis from a unit (bi)vector.
plotly:       true
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

\\
Toy C code can be found: [HERE](http://github.com/Marc-B-Reynolds/Stand-alone-junk/blob/master/src/Posts/ortho_basis.c)

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
  \mathbf{x_f}  & = & \left(  1+\frac{v_x^2}{v_z+1},~-\frac{v_x v_y}{v_z+1},~-v_x  \right) \label{frisvad_x} \\
  \mathbf{y_f}  & = & \left(  -\frac{v_x v_y}{v_z+1},~1-\frac{v_y^2}{v_z+1},~-v_y  \right) \label{frisvad_y} 
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
  float sz = -sgn(z);
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
  float sz = sgn(z);
  float a  = y/(fabsf(z)+1.f);
  float b  = y*a;
  float c  = x*a;
  
  vec3_set(xp, z+sz*b, sz*c, x); 
  vec3_set(yp, c,   1.f-b, -sz*y); 
}
{% endhighlight %}

<br>

------

Pixar's method <small>and inspiration thereof</small>
------

<small>addition added 20170401</small>

\\
A recently published paper from Pixar[^pixar] presents a new method directly derived from that of Frisvad $\eqref{frisvad_x} \eqref{frisvad_y}$ which is both efficient and high quality.  Translating their branch-free version into equations we have:

$$
\begin{eqnarray}
  \mathbf{x_p}  & = & 
\left(  1-\frac{s_z v_x^2}{v_z+s_z},~-\frac{s_z v_x v_y}{v_z+s_z},~-s_z v_x  \right) \label{pixar_x} \\
  \mathbf{y_p}  & = & 
\left(  -\frac{v_x v_y}{v_z+s_z},~s_z-\frac{v_y^2}{v_z+s_z},~-v_y  \right) \label{pixar_y} 
\end{eqnarray}
$$

\\
Over reals these are equivalent to $\eqref{xpn}$ and $\eqref{ypnn}$ (so it returns the opposite results of the *2a* method above). Their strategy for handling the negative half sphere is obviously superior to the one I used above. Where I was computing both abs and sign of z, they only use the latter.  The transform is simple: multiply top and bottom through by $s_z$:

$$
\frac{a}{\abs{v_z}+1} = \frac{s_z a}{v_z+s_z}
$$

\\
After the palmface I wanted to see what could be done with the new puzzle piece.

Starting with Pixar's method we can negate both the returned vectors and on X64 hardware this lowers the number of issues from 36 to 32[^gcc].  Additionally both GCC and clang fail to remove the redundant `sz*x` term which bring us to 30:

{% highlight c %}
void ortho_basis_pixar_r1(vec3_t* v, vec3_t* xp, vec3_t* yp)
{
  // -- snip start --
  float x  = v->x;
  float y  = v->y;
  float z  = v->z; 
  float sz = sgn(z);
  float a  = 1.0f/(sz+z);
  // -- snip end --
  float sx = sz*x;     // shouldn't be needed but is for gcc & clang, saves 2 issues
  float b  = x*y*a;
  vec3_set(xp, sx*x*a - 1.f, sz*b, sx);
  vec3_set(yp, b, y*y*a-sz, y);
}
{% endhighlight %}

\\
For X86 we've already beaten method 2a in terms of issues, which requires 32.

If we can allow left-handed results on one of the half-spheres then we can multiply the $\mathbf{x}'$ result of the previous through by $s_z$ which lowers the number of X64 issues to 28 (same as lower quality method 2 above): 

{% highlight c %}
void ortho_basis_pixar_l1(vec3_t* v, vec3_t* xp, vec3_t* yp)
{
  // -- cut-n-paste snipped out part --
  float b  = x*y*a;
  vec3_set(xp, x*x*a - sz, b, x);
  vec3_set(yp, b, y*y*a - sz, y);
}
{% endhighlight %}

\\
So far we haven't changed the error.  Since we are not allowing the compiler to treat FP as reals there are some redundant computations given the order of operations.  Choosing to pull out `y*a` drops the number of issues to 29/26 respectively and my RMS measures[^rms] show a minor negative impact[^stddiff] of ~0.003%. (Both compilers figure out the `sz*x` terms here)

{% highlight c %}
void ortho_basis_pixar_r2(vec3_t* v, vec3_t* xp, vec3_t* yp)
{
  // -- cut-n-paste snipped out part --
  float ya = y*a;
  float b  = x*ya;
  vec3_set(xp, sz*x*x*a - 1.f, sz*b, sz*x);
  vec3_set(yp, b, y*ya-sz, y);
}

void ortho_basis_pixar_l2(vec3_t* v, vec3_t* xp, vec3_t* yp)
{
  // -- cut-n-paste snipped out part --
  float ya = y*a;
  float b  = x*ya;
  vec3_set(xp, x*x*a - sz, b, x);
  vec3_set(yp, b, y*ya - sz, y);
}
{% endhighlight %}

\\
That exhasts the things that seem profitable and jump out at me to try working with x component reduction used in $\eqref{frisvad_x}$.

Returning the my previous choice of x-component expansion: take $\eqref{xpnn} \eqref{ypnn}$, change to Pixar's strategy and negating both give a left-handed system on one half sphere.

{% highlight c %}
void ortho_basis_l1(vec3_t* v, vec3_t* xp, vec3_t* yp)
{
  float x  = v->x;
  float y  = v->y;
  float z  = v->z; 
  float sz = sgn(z);
  float a  = y/(z+sz);
  float b  = y*a;
  float c  = x*a;
  
  vec3_set(xp, -z-b, c, x);
  vec3_set(yp, c, b-sz, y);
}
{% endhighlight %}

\\
and multiply $\mathbf{y}'$ by $s_z$ for right handed version:

{% highlight c %}
void ortho_basis_r1(vec3_t* v, vec3_t* xp, vec3_t* yp)
{
  float x  = v->x;
  float y  = v->y;
  float z  = v->z; 
  float sz = sgn(z);
  float a  = y/(z+sz);
  float b  = y*a;
  float c  = x*a;
  vec3_set(xp, -z-b, c, x);
  vec3_set(yp, sz*c, sz*b-1, sz*y);
}
{% endhighlight %}

\\
Note I'm only using X86 issues as a quick-n-dirty measure of goodness (easy to copy-n-past into godbolt).

A summary of results:


{: .center }
| method    | x86-64 |  RMS  | rel-RMS | notes |
|:----------|:------:|:-----:|:-------:|:-----:|
| basis_2   |     28 | 2.5868&times;10<sup>-8</sup> | -7.028% |  |
| basis_2a  |     32 | 2.5868&times;10<sup>-8</sup> | -7.028% |  |
| pixar     |     36 | 2.4170&times;10<sup>-8</sup> |    -   | as per paper |
| pixar_r1  |     30 | 2.4170&times;10<sup>-8</sup> |    -   |  |
| pixar_r2  |     29 | 2.4215&times;10<sup>-8</sup> | -0.003% |  |
| pixar_l1  |     28 | 2.4170&times;10<sup>-8</sup> |    -   |  |
| pixar_l2  |     26 | 2.4215&times;10<sup>-8</sup> | -0.003% |  |
| basis_r1  |     27 | 2.5868&times;10<sup>-8</sup> | -7.028% |  |
| basis_l2  |     24 | 2.5868&times;10<sup>-8</sup> | -7.028% |  |


<div id="error" style="width:100%"></div>

<br>


------

Other methods
------

\\
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
[^pixar]:    *"Building an Orthonormal Basis, Revisited"*, Tom Duff, James Burgess, Per Christensen, Christophe Hery, Andrew Kensler, Max Liani, and Ryusuke Villemin, Journal of Computer Graphics Techniques Vol. 6, No. 1, 2017 [JCGT](http://jcgt.org/published/0006/01/01/)
[^gcc]:      Instruction issue numbers are from x86-64 gcc 6.3 with -O3.  x86-64 clang 4.0.0 -O3 numbers appear to be the same.
[^rms]:      Take with a grain-of-salt. My RMS measures are not in agreement with those reported by Pixar.
[^stddiff]:  Computed the standard way:  $1-\frac{e}{e_0}$

<script>


var e1 = [8.346786e-08, 8.172310e-08, 8.259445e-08, 8.304380e-08, 8.365497e-08, 8.453741e-08, 8.428199e-08, 8.422571e-08, 8.570216e-08, 8.648852e-08, 8.651539e-08, 8.771895e-08, 8.783508e-08, 9.070049e-08, 9.612754e-08, 9.401791e-08, 9.505407e-08, 9.498568e-08, 9.460635e-08, 9.748003e-08, 9.872634e-08, 9.949827e-08, 9.799122e-08, 9.819922e-08, 1.025851e-07, 1.022257e-07, 1.029827e-07, 1.084380e-07, 1.073212e-07, 1.076929e-07, 1.071350e-07, 1.100539e-07, 1.092485e-07, 1.093655e-07, 1.110851e-07, 1.115556e-07, 1.137887e-07, 1.119455e-07, 1.165186e-07, 1.124847e-07, 1.233182e-07, 1.222818e-07, 1.219135e-07, 1.222025e-07, 1.241268e-07, 1.238678e-07, 1.209158e-07, 1.315551e-07, 1.247434e-07, 1.277628e-07, 1.312547e-07, 1.293994e-07, 1.327685e-07, 1.336069e-07, 1.334828e-07, 1.393003e-07, 1.349003e-07, 1.367431e-07, 1.386522e-07, 1.405115e-07, 1.431776e-07, 1.452278e-07, 1.445363e-07, 1.468288e-07, 1.512062e-07, 1.493344e-07, 1.461056e-07, 1.478839e-07, 1.543394e-07, 1.514666e-07, 1.545157e-07, 1.595809e-07, 1.577133e-07, 1.598156e-07, 1.627747e-07, 1.616733e-07, 1.641209e-07, 1.659786e-07, 1.705195e-07, 1.655401e-07, 1.695834e-07, 1.762140e-07, 1.814915e-07, 1.835591e-07, 1.846084e-07, 1.972617e-07, 1.866804e-07, 1.888271e-07, 1.881156e-07, 1.879487e-07, 1.917221e-07, 1.929102e-07, 1.941487e-07, 1.951785e-07, 2.022512e-07, 2.094967e-07, 2.037940e-07, 2.013958e-07, 2.054880e-07, 2.067102e-07, 2.067102e-07, 2.054880e-07, 2.013958e-07, 2.037940e-07, 2.094967e-07, 2.022512e-07, 1.951785e-07, 1.941487e-07, 1.929102e-07, 1.917221e-07, 1.879487e-07, 1.881156e-07, 1.888271e-07, 1.866804e-07, 1.972617e-07, 1.846084e-07, 1.835591e-07, 1.814915e-07, 1.762140e-07, 1.695834e-07, 1.655401e-07, 1.705195e-07, 1.659786e-07, 1.641209e-07, 1.616733e-07, 1.627747e-07, 1.598156e-07, 1.577133e-07, 1.595809e-07, 1.545157e-07, 1.514666e-07, 1.543394e-07, 1.478839e-07, 1.461056e-07, 1.493344e-07, 1.512062e-07, 1.468288e-07, 1.445363e-07, 1.452278e-07, 1.431776e-07, 1.405115e-07, 1.386522e-07, 1.367431e-07, 1.349003e-07, 1.393003e-07, 1.334828e-07, 1.336069e-07, 1.327685e-07, 1.293994e-07, 1.312547e-07, 1.277628e-07, 1.247434e-07, 1.315551e-07, 1.209158e-07, 1.238678e-07, 1.241268e-07, 1.222025e-07, 1.219135e-07, 1.222818e-07, 1.233182e-07, 1.124847e-07, 1.165186e-07, 1.119455e-07, 1.137887e-07, 1.115556e-07, 1.110851e-07, 1.093655e-07, 1.092485e-07, 1.100539e-07, 1.071350e-07, 1.076929e-07, 1.073212e-07, 1.084380e-07, 1.029827e-07, 1.022257e-07, 1.025851e-07, 9.819922e-08, 9.799122e-08, 9.949827e-08, 9.872634e-08, 9.748003e-08, 9.460635e-08, 9.498568e-08, 9.505407e-08, 9.401791e-08, 9.612754e-08, 9.070049e-08, 8.783508e-08, 8.771895e-08, 8.651539e-08, 8.648852e-08, 8.570216e-08, 8.422571e-08, 8.428199e-08, 8.453741e-08, 8.365497e-08, 8.304380e-08, 8.259445e-08, 8.172310e-08, 8.346786e-08];

var e2 = [5.421286e-08, 5.529221e-08, 5.492030e-08, 5.668010e-08, 5.581707e-08, 5.656291e-08, 5.655056e-08, 5.652862e-08, 5.755258e-08, 5.827501e-08, 5.825030e-08, 6.026310e-08, 5.949323e-08, 6.121726e-08, 6.501582e-08, 6.282821e-08, 6.539093e-08, 6.521670e-08, 6.828747e-08, 6.831470e-08, 7.004165e-08, 6.864687e-08, 7.016793e-08, 7.135286e-08, 7.189237e-08, 7.351777e-08, 7.522851e-08, 7.621363e-08, 7.782879e-08, 8.058189e-08, 8.205831e-08, 7.986833e-08, 8.073583e-08, 8.699194e-08, 8.565545e-08, 8.379074e-08, 8.768495e-08, 8.836715e-08, 8.940996e-08, 8.559652e-08, 9.045983e-08, 9.253621e-08, 9.128107e-08, 9.113682e-08, 9.151952e-08, 9.445538e-08, 9.727483e-08, 9.465918e-08, 9.449763e-08, 9.661785e-08, 9.462272e-08, 9.820104e-08, 9.727161e-08, 9.683593e-08, 9.737888e-08, 9.714789e-08, 9.932323e-08, 9.853940e-08, 1.021374e-07, 1.017948e-07, 1.017019e-07, 1.053757e-07, 1.062081e-07, 1.072450e-07, 1.067818e-07, 1.087370e-07, 1.106508e-07, 1.144913e-07, 1.124766e-07, 1.158639e-07, 1.131636e-07, 1.169442e-07, 1.177146e-07, 1.219969e-07, 1.213930e-07, 1.246513e-07, 1.252047e-07, 1.262536e-07, 1.255877e-07, 1.259231e-07, 1.270965e-07, 1.287834e-07, 1.395258e-07, 1.336800e-07, 1.340486e-07, 1.428696e-07, 1.387385e-07, 1.458600e-07, 1.441353e-07, 1.435760e-07, 1.442086e-07, 1.514289e-07, 1.474636e-07, 1.549530e-07, 1.542205e-07, 1.548432e-07, 1.589722e-07, 1.536757e-07, 1.601940e-07, 1.633581e-07, 1.633581e-07, 1.601940e-07, 1.536757e-07, 1.558218e-07, 1.548432e-07, 1.542205e-07, 1.549530e-07, 1.474636e-07, 1.514289e-07, 1.442086e-07, 1.435760e-07, 1.441353e-07, 1.342940e-07, 1.387385e-07, 1.388419e-07, 1.340486e-07, 1.336800e-07, 1.395258e-07, 1.287834e-07, 1.270965e-07, 1.259231e-07, 1.255877e-07, 1.262536e-07, 1.252047e-07, 1.246513e-07, 1.213930e-07, 1.219969e-07, 1.177146e-07, 1.169442e-07, 1.131636e-07, 1.158639e-07, 1.124766e-07, 1.144913e-07, 1.106508e-07, 1.087370e-07, 1.067818e-07, 1.072450e-07, 1.062081e-07, 1.053757e-07, 1.017019e-07, 1.017948e-07, 1.021374e-07, 9.853940e-08, 9.932323e-08, 9.714789e-08, 9.737888e-08, 9.683593e-08, 9.727161e-08, 9.820104e-08, 9.462272e-08, 9.661785e-08, 9.449763e-08, 9.465918e-08, 9.727483e-08, 9.445538e-08, 9.151952e-08, 9.113682e-08, 9.128107e-08, 9.253621e-08, 9.045983e-08, 8.559652e-08, 8.940996e-08, 8.836715e-08, 8.768495e-08, 8.379074e-08, 8.565545e-08, 8.699194e-08, 8.073583e-08, 7.986833e-08, 8.205831e-08, 8.058189e-08, 7.782879e-08, 7.621363e-08, 7.522851e-08, 7.351777e-08, 7.189237e-08, 7.135286e-08, 7.016793e-08, 6.864687e-08, 7.004165e-08, 6.831470e-08, 6.828747e-08, 6.521670e-08, 6.539093e-08, 6.282821e-08, 6.501582e-08, 6.121726e-08, 5.949323e-08, 6.026310e-08, 5.825030e-08, 5.827501e-08, 5.755258e-08, 5.652862e-08, 5.655056e-08, 5.656291e-08, 5.581707e-08, 5.668010e-08, 5.492030e-08, 5.529221e-08, 5.421286e-08];

var e3 = [5.525352e-08, 5.503035e-08, 5.521761e-08, 5.574448e-08, 5.665871e-08, 5.726875e-08, 5.692867e-08, 5.754394e-08, 5.747906e-08, 5.814501e-08, 5.858978e-08, 5.957587e-08, 5.915469e-08, 6.175146e-08, 6.308243e-08, 6.354516e-08, 6.362824e-08, 6.755210e-08, 6.527229e-08, 6.709297e-08, 6.617284e-08, 6.981048e-08, 6.913312e-08, 7.295046e-08, 7.509046e-08, 7.302444e-08, 7.568009e-08, 7.577605e-08, 7.917795e-08, 8.193540e-08, 8.344853e-08, 8.234191e-08, 8.307613e-08, 8.416659e-08, 8.399125e-08, 8.512576e-08, 8.494308e-08, 8.840386e-08, 8.713028e-08, 9.154841e-08, 8.765430e-08, 9.041441e-08, 9.025087e-08, 9.222824e-08, 9.338154e-08, 9.566909e-08, 9.549178e-08, 9.504764e-08, 9.778892e-08, 9.849469e-08, 9.570064e-08, 1.010246e-07, 9.619046e-08, 9.698321e-08, 1.019651e-07, 9.730119e-08, 1.025287e-07, 1.012587e-07, 9.982812e-08, 1.007528e-07, 1.066140e-07, 1.054583e-07, 1.073725e-07, 1.114864e-07, 1.146415e-07, 1.110054e-07, 1.123303e-07, 1.131016e-07, 1.191513e-07, 1.156238e-07, 1.229841e-07, 1.226790e-07, 1.207866e-07, 1.223726e-07, 1.217327e-07, 1.247132e-07, 1.262953e-07, 1.269429e-07, 1.351482e-07, 1.280721e-07, 1.304202e-07, 1.324718e-07, 1.422551e-07, 1.338617e-07, 1.393933e-07, 1.359860e-07, 1.440760e-07, 1.394110e-07, 1.420514e-07, 1.400995e-07, 1.535680e-07, 1.465763e-07, 1.499561e-07, 1.556371e-07, 1.587793e-07, 1.533983e-07, 1.585293e-07, 1.562631e-07, 1.555635e-07, 1.576881e-07, 1.576881e-07, 1.555635e-07, 1.544795e-07, 1.585293e-07, 1.533983e-07, 1.521094e-07, 1.530449e-07, 1.499561e-07, 1.465763e-07, 1.441623e-07, 1.400995e-07, 1.420514e-07, 1.394110e-07, 1.440760e-07, 1.359860e-07, 1.393933e-07, 1.338617e-07, 1.422551e-07, 1.324718e-07, 1.304202e-07, 1.277102e-07, 1.351482e-07, 1.269429e-07, 1.262953e-07, 1.247132e-07, 1.217327e-07, 1.223726e-07, 1.207866e-07, 1.226790e-07, 1.229841e-07, 1.156238e-07, 1.191513e-07, 1.131016e-07, 1.123303e-07, 1.110054e-07, 1.146415e-07, 1.114864e-07, 1.073725e-07, 1.054583e-07, 1.066140e-07, 1.007528e-07, 9.982812e-08, 1.012587e-07, 1.025287e-07, 9.730119e-08, 1.019651e-07, 9.698321e-08, 9.619046e-08, 1.010246e-07, 9.570064e-08, 9.849469e-08, 9.778892e-08, 9.504764e-08, 9.549178e-08, 9.566909e-08, 9.338154e-08, 9.222824e-08, 9.025087e-08, 9.041441e-08, 8.765430e-08, 9.154841e-08, 8.713028e-08, 8.840386e-08, 8.494308e-08, 8.512576e-08, 8.399125e-08, 8.416659e-08, 8.307613e-08, 8.234191e-08, 8.344853e-08, 8.193540e-08, 7.917795e-08, 7.577605e-08, 7.568009e-08, 7.302444e-08, 7.509046e-08, 7.295046e-08, 6.913312e-08, 6.981048e-08, 6.617284e-08, 6.709297e-08, 6.527229e-08, 6.755210e-08, 6.362824e-08, 6.354516e-08, 6.308243e-08, 6.175146e-08, 5.915469e-08, 5.957587e-08, 5.858978e-08, 5.814501e-08, 5.747906e-08, 5.754394e-08, 5.692867e-08, 5.726875e-08, 5.665871e-08, 5.574448e-08, 5.521761e-08, 5.503035e-08, 5.525352e-08];

var etrace0 = { y: e1, x0:-1, dx:1/100, name: 'bar'};
var etrace1 = { y: e2, x0:-1, dx:1/100, name: 'foo'};
var etrace2 = { y: e3, x0:-1, dx:1/100, name: 'baz'};

var edata   = [etrace0,etrace1,etrace2];

var elayout = {
  xaxis: {title: "RMS" }
};

Plotly.newPlot('error', edata, elayout);

</script>
