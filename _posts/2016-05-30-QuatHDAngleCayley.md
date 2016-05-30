---
layout:       post
title:        Quaternion half/double angle and Cayley transforms
categories:   [quaternions]
tags:         [interpolation, quantization, domain reduction]
description:  a note on two transform pairs which reduce the working domain in quaternion space
plotly:       true
---

<div class="alert alert-warning" role="alert" markdown="1">
DRAFT: come back later
</div>

This post is about two transform pairs of quaternion that reduce the working domain which are not commonly discussed.  Specifically they both reduce the implied angle in quaternion space.  This allows storing a quaternion value with three elements more attractive, but more importantly reduces the range of the implied angle. This angle reduction can be useful for both interpolation and quantization for example.  The usefulness of these maps depend on the target error bound and maximum implied angle the system allows.  There are relatively well known techniques which use unit magnitude and that $Q$ and $-Q$ are equivalent for alternate storage and/or quantization.  I will later examine the case of interpolation and will focus on the basic properties and present some simple quantization.

Once we get past the preliminary junk we'll be almost entirely talking about a 2D plane and complex numbers.  I will also drop most of the formalism.  I'm also going to skip on actually defining the sets and will simply note how I name them and they are hopefully obvious by context.

I'm using pseudo-GLSL for code snippets in an attempt at clarity.  I will make very few nods to optimizations other than carrying through some of the derivations.

------

Preliminary stuff
------

\\
We generally work with quaternions $\left(\mathbb{H}\right)$ represented by four real values:

$$ \begin{equation} \label{qxform}
Q = a + \left(x,~y,~z\right) 
\end{equation} $$

\\
Where I will call "$a$" the scalar part and $(x,y,z)$ the bivector part. For our purposes here we can consider a 3D bivector to be equivalent to a 3D vector.

We will also want to express a quaternion in terms of some implied unit bivector $\mathbf{u}$ $\left( \mathbb{S}^{2} \right) $ with implied magnitude $b$:

$$ \begin{equation} \label{eq:qsb}
Q  =  a + b~\mathbf{u}
\end{equation} $$ 

\\
and in polar form[^polar] with an implied angle (measurement in quaternion space) of $\theta$:

$$ \begin{equation} \label{eq:qpolar}
Q = \cos\left(\theta\right)+\sin\left(\theta\right)~\mathbf{u}
\end{equation} $$ 

\\
Since all non-zero scalar multiples of a quaternion represent the same rotation (projective line $sQ$, for unit quaternion $Q$ and non-zero scalar $s$) so we can restrict ourselves to unit quaternions.  Equation $\eqref{eq:qpolar}$ already makes this assumption.  This gives the four dimension unit sphere $\left( \mathbb{S}^{3} \right) $ as a working domain.

The projective line intersects the unit-sphere at $Q$ and $-Q$, so we can further restrict ourself to quaternions with scalars part that are zero or positive $\left( a \geq 0 \right)$.  In axis-angle think this is negating the angle and reversing the axis.  The reduces the working domain to the half sphere $\left( \mathbb{S}^{3+} \right) $. $ \newcommand{\hcomplex}[1]{\mathbb{C}_{#1}} $  So far all standard fare.

If we take the set of all quaternion with fixed $\mathbf{u}$ definable in forms $\eqref{eq:qsb}~\eqref{eq:qpolar}$ they form a plane which spans the line-of-reals and $\mathbf{u}$.  [This plane forms a commutative subalgebra of $\mathbb{H}$ isomorphic to $\mathbb{C}$]({{site.base}}/quaternions/2016/05/17/QuatAsComplex.html#complexPlane), we can denote this plane as $\hcomplex{u}$ or as the {1,Q}-plane.  The set of all quaternions is the collection of all unique complex planes and they intersect at the line of reals.

The transform pairs which we will be considering are what I call [*(quaternion valued) complex functions*]({{site.base}}/quaternions/2016/05/17/QuatAsComplex.html#complexFunction) which simply means if we have a function $f$ which given input in some complex plane then the result is in the same complex plane.  In other words we are working in two dimensions using the standard set of complex numbers.

We are at the following in terms of $\eqref{eq:qsb}$:

* effectively working with the complex number: $z = a + b~\mathbf{i} $
* unit magnitude: $ a^2 + b^2 =1 $
* positive scalar: $ a $ on $\left[0,1\right] $
* $ b $ on $\left[0,1\right] $, since $b=\sqrt{x^2+y^2+z^2}$. Any sign of $b$ is contained by the implied $\mathbf{u}$
* the implied angle $\theta$ is therefore on $\left[0,~\frac{\pi}{2}\right]$
* the domain is the arc of unit circle in the first quadrant, axes included.

The presented forward transforms are also in the first quadrant, axes included.

<br>

------

Half/double angle
------
{:#halfDoubleAngle}

\\
Quaternions naturally contain a half-angle/double-angle relation.  To convert from an axis-angle representation we halve the angle of the rotation.  If we have a rotation represented in axis-angle form with the axis represent by the unit vector $\hat{u}$ and an angle of rotation $\alpha$, then the conversion to quaternions can be expressed as:

$$ \begin{equation*}
\cos\left(\frac{\alpha}{2}\right)+\sin\left(\frac{\alpha}{2}\right)~\mathbf{u}
\end{equation*} $$ 

\\
where the unit vector $\hat{u}$ and unit bivector $u$ are the same $(x,y,z)$ values.  In turn this half-angle representation is doubled when we apply the similarity transform[^rotViz] (the rotation matrix conversion is the similarity transform applied to an orthonormal basis).

$$ \begin{eqnarray*}
A & = & QAQ^{-1} \\
  & = & QAQ^{*} 
\end{eqnarray*} $$


We can introduce a second half-angle/double-angle transform pair.  In complex number these translate to square root and squaring (or trig identities or geometry).  The half angle transform is:

$$ \begin{eqnarray*}
  \sqrt{z} & = & \sqrt{a + b~\mathbf{i}} \\
           & = & \frac{1}{\sqrt{2+2a}} \left(1+a + b~\mathbf{i} \right)
\end{eqnarray*} $$

\\
And the double-angle transform:

$$ \begin{eqnarray}
  z^2 & = & \left(a + b~\mathbf{i}\right)^2 \nonumber \\
      & = & a^2 - b^2 + 2ab~\mathbf{i}      \nonumber \\
      & = & 1-2a^2 + 2ab~\mathbf{i}         \label{zsq}  \\
      & = & 1-2b^2 + 2b\sqrt{1-b^2}~\mathbf{i}         \label{zsqBv}
\end{eqnarray} $$


\\
Equation $\eqref{zsq}$ is for reconstructing from all four components and $\eqref{zsqBv}$ for when we only have the values of the bivector part.

After the half-angle transform we have:

* $ a^2 + b^2 =1 $
* $\theta \in \left[0,~\frac{\pi}{4}\right]$
* $ a \in \left[\frac{1}{\sqrt{2}},1\right] $
* $ b \in \left[0, \frac{1}{\sqrt{2}}\right] $
* $ a \ge b $

\\
For the complex plane as a whole we are halving the angle with respect to the line of reals and the magnitude is square-rooted.  This maps the negative half-plane to the positive.  Steven Wittens' blog post ["How to Fold a Julia Fractal"](http://acko.net/blog/how-to-fold-a-julia-fractal/) has some handy [mathbox](http://gitgud.io/unconed/mathbox) visualizations in  "Like Hands on a Clock" section.  From step 29->30 shows the square-root (with negative $x$ pointing down).

Remember in the plane all points along the projective line are equivalent.  We can think of the projective line as being two half-lines on either size of the origin that have an angle of $\pi$ between them.  After the half-angle transform these half lines have a relative angle of $\frac{\pi}{2}$ and since nothing has happened to the properties in the plane, so we now have two orthogonal projective lines all of which are equivalent representations.  That means four unit quaternions are equivalent. ("Like Hands on a Clock": step 30->31).  When we square again these two images of the projective line become one again.

The half-angle transform does not effect angle measurements between bivectors so we cannot simply treat half-angle representations as the standard and expect everything to work.  As an example we cannot simply compose rotations with a product. In this additional half-angle representation then the proper rotation is $ AABB $.  If we attempted to compose with product and then square we would have $ ABAB $.  I will say more about this in later posts.

<br>

{% highlight c %}
// quaternion 'q' as vec4: q.w = scalar part 'a' , q.xyz = bivector part 'b'

// sqrt(a + bu)
vec4 q_usqrt(vec4 q) {
  float d = 1.0+q.w;           // 1+a
  float s = inversesqrt(d+d);  // 1/sqrt(2+2a)
  vec3  b = s*q.xyz;           // b/sqrt(2+2a) u
  return vec4(b, d*s);         // (1+a)/sqrt(2+2a) + b/sqrt(2+2a) u
}

// (a + bu)^2
vec4 q_upow2(vec4 q) {
  float s = q.w+q.w;           // 2a
  float w = s*q.w;             // 2a^2
  vec3  b = s*q.xyz;           // 2ab u
  return vec4(b, 1.0-w);       // (1-2a^2) + 2ab u
}

// (a + bu)^2  where we only have the 'bu' term
vec4 q_upow2b(vec3 v) {
  float d = dot(v,v);               // b^2
  float s = sqrt(1.0-d);            // sqrt(1-b^2)
  return vec4((s+s)*v, 1.0-(d+d));  // (1-2b^2) + 2b sqrt(1-b^2) u
}

{% endhighlight %}

<br>

------

Log/exp transform
------
{:#expLog}

\\
Another pair of complex functions are $\log$ and $\exp$.  For unit quaterions (complex numbers) $\log$ can be expressed in terms of the polar form \eqref{eq:qpolar} as:

$$ \begin{eqnarray*}
  \log\left( \cos\left(\theta \right) + \sin\left(\theta \right)\mathbf{i}   \right) & = &  \theta~\mathbf{i} \\
\end{eqnarray*} $$

\\
which maps a unit quaternion to a bivector.  The reverse transform $\exp$ when given a bivector can be expressed as:

$$ \begin{eqnarray*}
  \exp\left(\theta~\mathbf{i}\right) & = & \cos\left(\theta \right) + \sin\left(\theta \right)\mathbf{i}
\end{eqnarray*} $$

\\
which is simply another (and common) way to express the polar form.  Notice that $\log$ is implicitly the conversion from an axis-angle representation and is sometime called *Euler form*.  Basically it is simply a connection between the algebra and trig operations.  Since we are on a limited range for both transforms we could use a combination of trig identities and function approximation to come up with reasonably accurate and fast versions, however here we are only going to consider these as the "ideal target" transforms since $\log$ perfectly linearizes the angle where we have:

* $ a = 0 $
* $ b \in \left[0, \frac{\pi}{2}\right] $

\\
For the space as a whole we are mapping our 4D unit half-sphere input to a 3D ball with a radius of $\frac{\pi}{2}$.

<br>

------

Cayley transform
------
{:#cayley}

\\
The Cayley transform dates back to the 1840s.  It is a special case of a stereographic projection and it is a Mobius transformation (linear fractional transformation).  It is also the connection between the Lie algebra and group and is related to tangent half-angle substitution.  If we call the forward transform $f$ and the reverse $g$ then we have:

$$ \begin{eqnarray*}
f\left(z\right) & = &  \frac{z-1}{z+1} \\
                & = &  \left(z-1\right)\left(z+1\right)^{-1}                 \\
                & = &  \frac{\left(z-1\right)\left(z^*+1\right)}{\left(z+1\right)\cdot\left(z+1\right)}       \\
		\\
g\left(z\right) & = &  \frac{1+z}{1-z} = \frac{\left(1+z\right)\left(1-z^*\right)}{\left(1-z\right)\cdot\left(1-z\right)}
\end{eqnarray*} $$

\\
Expanding the forward transform with $\eqref{eq:qsb}$:

$$ \begin{eqnarray*}
f\left(a+b~\mathbf{i}\right) & = &
                \frac{\left(a-1+b~\mathbf{i}\right)\left(a+1-b~\mathbf{i}\right)}{\left(a+1+b~\mathbf{i}\right)\cdot\left(a+1+b~\mathbf{i}\right)} \\
		& = & \frac{2b}{\left(a+1\right)^2+b^2}~\mathbf{i} \\
                & = & \frac{2b}{2\left(1+a\right)}~\mathbf{i} \\
		& = & \frac{b}{1+a}~\mathbf{i}
\end{eqnarray*} $$

\\
Like $\log$ the result is a bivector.  In polar form we have:

$$
f\left(\cos\left(\theta\right)+\sin\left(\theta\right)~\mathbf{i}\right)  =  \frac{\sin\left(\theta\right)}{\cos\left(\theta\right)+1}~\mathbf{i}
 =  \tan\left(\frac{\theta}{2}\right)~\mathbf{i}
$$

\\
We only need the reverse transform for bivector input:

$$ \begin{eqnarray*}
g\left(b~\mathbf{i}\right) & = &  \frac{\left(1+b~\mathbf{i}\right)\left(1+b~\mathbf{i}\right)}{\left(1-b~\mathbf{i}\right)\cdot\left(1-b~\mathbf{i}\right)} \\
                & = &  \frac{1-b^2+2b~\mathbf{i}}{1+b^2} \\
                & = &  \frac{2}{1+b^2} - 1 + \frac{2}{1+b^2}b~\mathbf{i}
\end{eqnarray*} $$

\\
After the forward transform we have:

* $ a = 0 $
* $ b \in \left[0, 1\right] $

For the space as a whole we are mapping our 4D unit half-sphere input to the 3D unit ball.  If we were to apply the forward transform above to the negative 4D half-sphere it maps to all of the 3D space outside of the unit ball.  Inverting the original transforms is the opposite transform pair (negative to unit ball, positive to outside). So we have two equivalent coordinates in this space as well $Q$ and $Q^{-1}$.

<br>

{% highlight c %}

vec3 q_cayley(vec4 q) { return q.xyz*(1.0/(1.0+q.w)); }

vec4 q_cayley(vec3 v) {
  float s = 2.0/(1.0+dot(v,v)); // 2/(1+b^2)
  return vec4(s*v, s-1.0);
}

{% endhighlight %}

<br>

------

Bivector part
------
{:#bivectorPart}

None of these transforms effect the implied unit bivector $u$. XXX
<br>

$$ \begin{eqnarray*}
 x \in \left[\frac{1}{\sqrt{3}},~1\right] & \approx & \left[0.57735,~1\right] \\
 y \in \left[0,\frac{1}{\sqrt{2}}\right]  & \approx & \left[0,~.70711\right] \\
 z \in \left[0,\frac{1}{\sqrt{3}}\right] & \approx & \left[0,~0.57735\right]
\end{eqnarray*} $$

XXX


<br>

------

Comparison of mappings
------
{:#compare}

\\
These figures use [plotly.js](http://plot.ly/javascript/). You can click on the function in the legend to hide/show that plot.

\\
The first figure is a simply comparison of how the magnitude of the bivector behaves with respect to the implied angle.  The $x$-axis is the normalized angle, so in terms of the represented angle of rotation in 3D we would multiply by $\pi$.  The $y$-axis is normalized magnitude of bivector (maximum value of $b$ divided out).  If we were to put an upper limit on the implied angle then these plots would grow closer to that of $\log$.

<div id="fig1" style="width:100%"></div>

\\
The second figure is the same as the first with simply the linear part subtracted out.  If we were to put an upper limit on the implied angle then the shape of plots stay pretty much the same.

<div id="fig2" style="width:100%"></div>

\\
We can think of the first figure as being plots of cumulative density where the $y$ value is "how many" and the corresponding $x$ is "where".  As an example if we are linearly quantizing ($N$ samples) and examine $x$ at 0.5 then the $y$ value of each plot is a multiplier $m$ and $mN$ is the number of samples on normalized angle range [0, 0.5].  Or examine the $y$ at 0.5 and then $x$ tells us how far along the angle range the first half of samples represents.

\\
By taking the derivative we can examine localized density which gives a much better picture of how samples are distributed across the full range.  Like the second figure the shape of the plots do not change much when adding a cap to maximum implied angle.

<div id="fig3" style="width:100%"></div>

<br>

------

Strawman examples
------


<br>

------

References and Footnotes
------

[^polar]: Some text prefer $e^{u\theta} $. 
[^rotViz]: Quaternion rotation visualization: ([local post]({{site.base}}/quaternions/2016/05/17/VisualizeQuatRot.html))
[^plotly]: Some text prefer $e^{u\theta} $. 



<script>
var plot0 = {
x:[2.04082e-8, 0.000306718, 0.000613415, 0.00122681, 0.0024536, 
0.00490718, 0.00981434, 0.0196287, 0.0409084, 0.0607779, 0.0802576, 
0.101388, 0.121109, 0.142481, 0.163463, 0.183035, 0.204257, 0.22407, 
0.243493, 0.264567, 0.284231, 0.305546, 0.326471, 0.345985, 0.367151, 
0.386907, 0.408314, 0.429331, 0.448938, 0.470196, 0.490044, 0.509502, 
0.530611, 0.55031, 0.57166, 0.59262, 0.61217, 0.633371, 0.653162, 
0.672563, 0.693616, 0.713258, 0.734551, 0.754434, 0.773927, 0.795071, 
0.814805, 0.83619, 0.857186, 0.876771, 0.877103, 0.877435, 0.878098, 
0.879426, 0.88208, 0.887389, 0.898007, 0.898317, 0.898627, 0.899246, 
0.900486, 0.902964, 0.90792, 0.90823, 0.90854, 0.909159, 0.910399, 
0.912877, 0.917833, 0.918137, 0.918441, 0.919048, 0.920263, 0.922692, 
0.927552, 0.927855, 0.928159, 0.928766, 0.929981, 0.932411, 0.93727, 
0.937599, 0.937929, 0.938588, 0.939906, 0.942542, 0.947813, 0.948143, 
0.948472, 0.949131, 0.950449, 0.953085, 0.958357, 0.958665, 0.958972, 
0.959587, 0.960817, 0.963276, 0.963584, 0.963891, 0.964506, 0.965736, 
0.968196, 0.968503, 0.968811, 0.969425, 0.970655, 0.973115, 0.973422, 
0.97373, 0.974345, 0.975575, 0.978034, 0.978377, 0.978721, 0.979407, 
0.98078, 0.981123, 0.981466, 0.982153, 0.983526, 0.983869, 0.984212, 
0.984899, 0.986271, 0.986615, 0.986958, 0.987644, 0.989017, 0.98936, 
0.989704, 0.99039, 0.990733, 0.991076, 0.991763, 0.992106, 0.992449, 
0.993136, 0.993479, 0.993822, 0.994509, 0.994852, 0.995195, 0.995538, 
0.995881, 0.996225, 0.996568, 0.996911, 0.997254, 0.997597, 0.997941, 
0.998284, 0.998627, 0.99897, 0.999314, 0.999657, 1.],
y:[3.20571e-8, 0.000481791, 0.000963551, 0.00192707, 0.0038541, 
0.0077081, 0.0154157, 0.0308277, 0.0642145, 0.0953247, 0.125735, 
0.158588, 0.189092, 0.221945, 0.253955, 0.283565, 0.31537, 0.344746, 
0.37322, 0.403721, 0.431783, 0.461735, 0.490635, 0.517111, 0.545278, 
0.571024, 0.5983, 0.624422, 0.648179, 0.67324, 0.695962, 0.717582, 
0.740277, 0.760722, 0.782058, 0.802149, 0.820105, 0.838703, 0.855225, 
0.87062, 0.886409, 0.900267, 0.914322, 0.926522, 0.937607, 0.948636, 
0.957985, 0.967078, 0.974943, 0.981324, 0.981424, 0.981524, 0.981723, 
0.982118, 0.982894, 0.984396, 0.987194, 0.987271, 0.987349, 0.987502, 
0.987807, 0.988406, 0.989558, 0.989628, 0.989698, 0.989837, 0.990112, 
0.99065, 0.991682, 0.991744, 0.991805, 0.991926, 0.992166, 0.992636, 
0.993532, 0.993586, 0.993639, 0.993746, 0.993958, 0.994369, 0.995149, 
0.9952, 0.99525, 0.995351, 0.995548, 0.99593, 0.996642, 0.996684, 
0.996726, 0.996809, 0.996972, 0.997286, 0.997861, 0.997893, 0.997924, 
0.997986, 0.998106, 0.998337, 0.998364, 0.998392, 0.998446, 0.998552, 
0.998752, 0.998776, 0.9988, 0.998847, 0.998938, 0.999108, 0.999129, 
0.999149, 0.999188, 0.999264, 0.999405, 0.999423, 0.999441, 0.999477, 
0.999544, 0.99956, 0.999576, 0.999607, 0.999665, 0.999679, 0.999693, 
0.999719, 0.999767, 0.999779, 0.99979, 0.999812, 0.999851, 0.99986, 
0.999869, 0.999886, 0.999894, 0.999902, 0.999916, 0.999923, 0.99993, 
0.999942, 0.999948, 0.999953, 0.999963, 0.999967, 0.999972, 0.999975, 
0.999979, 0.999982, 0.999985, 0.999988, 0.999991, 0.999993, 0.999995, 
0.999996, 0.999998, 0.999999, 0.999999, 1., 1.],
mode: 'lines',
name: 'standard'
};

var plot1 = {
x:[2.04082e-8, 0.000306718, 0.000613415, 0.00122681, 0.0024536, 
0.00490718, 0.00981434, 0.0196287, 0.0409084, 0.0607779, 0.0802576, 
0.101388, 0.121109, 0.142481, 0.163463, 0.183035, 0.204257, 0.22407, 
0.243493, 0.264567, 0.284231, 0.305546, 0.326471, 0.345985, 0.367151, 
0.386907, 0.408314, 0.429331, 0.448938, 0.470196, 0.490044, 0.509502, 
0.530611, 0.55031, 0.57166, 0.59262, 0.61217, 0.633371, 0.653162, 
0.672563, 0.693616, 0.713258, 0.734551, 0.754434, 0.773927, 0.795071, 
0.814805, 0.83619, 0.857186, 0.876771, 0.898007, 0.917833, 0.93727, 
0.958357, 0.978034, 0.978377, 0.978721, 0.979407, 0.98078, 0.983526, 
0.989017, 0.98936, 0.989704, 0.99039, 0.991763, 0.994509, 0.994852, 
0.995195, 0.995881, 0.997254, 0.997597, 0.997941, 0.998627, 0.99897, 
0.999314, 0.999657, 1.],
y:[2.26678e-8, 0.000340678, 0.000681333, 0.00136264, 0.00272526, 
0.00545049, 0.0109009, 0.0218011, 0.0454299, 0.0674816, 0.0890848, 
0.112495, 0.134316, 0.157926, 0.181063, 0.202601, 0.225901, 0.247597, 
0.268807, 0.29175, 0.313085, 0.336128, 0.358657, 0.379581, 0.402175, 
0.423163, 0.445789, 0.467882, 0.488377, 0.510467, 0.530963, 0.550932, 
0.572449, 0.592387, 0.613836, 0.634725, 0.654054, 0.674841, 0.694077, 
0.712771, 0.732869, 0.75144, 0.77137, 0.789785, 0.807652, 0.826819, 
0.844501, 0.863434, 0.881785, 0.898687, 0.916774, 0.93343, 0.949538, 
0.966765, 0.9826, 0.982874, 0.983148, 0.983696, 0.984791, 0.986978, 
0.991337, 0.991609, 0.991881, 0.992424, 0.99351, 0.995678, 0.995948, 
0.996219, 0.99676, 0.997841, 0.998111, 0.998381, 0.998921, 0.999191, 
0.999461, 0.99973, 1.],
mode: 'lines',
name: 'half angle'
};

var plot2 = {
x:[2.04082e-8, 0.000306718, 0.000613415, 0.00122681, 0.0024536, 
0.00490718, 0.00981434, 0.0196287, 0.0409084, 0.0607779, 0.0802576, 
0.101388, 0.121109, 0.142481, 0.163463, 0.183035, 0.204257, 0.22407, 
0.243493, 0.264567, 0.284231, 0.305546, 0.326471, 0.345985, 0.367151, 
0.386907, 0.408314, 0.429331, 0.448938, 0.470196, 0.490044, 0.509502, 
0.530611, 0.55031, 0.57166, 0.59262, 0.61217, 0.633371, 0.653162, 
0.672563, 0.693616, 0.713258, 0.734551, 0.754434, 0.773927, 0.795071, 
0.814805, 0.83619, 0.857186, 0.876771, 0.898007, 0.917833, 0.93727, 
0.958357, 0.978034, 0.978377, 0.978721, 0.979407, 0.98078, 0.983526, 
0.989017, 0.98936, 0.989704, 0.99039, 0.991763, 0.994509, 0.994852, 
0.995195, 0.995881, 0.997254, 0.997597, 0.997941, 0.998627, 0.99897, 
0.999314, 0.999657, 1.],
y:[1.60285e-8, 0.000240896, 0.000481775, 0.000963535, 0.00192706, 
0.00385411, 0.00770832, 0.0154175, 0.0321404, 0.0477711, 0.0631178, 
0.0797991, 0.0954068, 0.112374, 0.129093, 0.144754, 0.161814, 
0.177824, 0.193605, 0.210833, 0.227018, 0.24469, 0.26218, 0.278628, 
0.296628, 0.313589, 0.332154, 0.350585, 0.367972, 0.387048, 0.405082, 
0.422984, 0.442669, 0.461302, 0.481799, 0.502246, 0.521625, 0.542994, 
0.563294, 0.583542, 0.605925, 0.627216, 0.650768, 0.673227, 0.69571, 
0.720645, 0.744461, 0.770898, 0.797529, 0.823009, 0.851378, 0.878599, 
0.906018, 0.936637, 0.966078, 0.966599, 0.967121, 0.968165, 0.970256, 
0.974451, 0.982895, 0.983425, 0.983956, 0.985017, 0.987144, 0.991411, 
0.991946, 0.992481, 0.993551, 0.995696, 0.996233, 0.99677, 0.997846, 
0.998384, 0.998922, 0.999461, 1.],
mode: 'lines',
name: 'Cayley'
};

var plot3 = {
x:[2.04082e-8, 0.000306718, 0.000613415, 0.00122681, 0.0024536, 
0.00490718, 0.00981434, 0.0196287, 0.0409084, 0.0607779, 0.0802576, 
0.101388, 0.121109, 0.142481, 0.163463, 0.183035, 0.204257, 0.22407, 
0.243493, 0.264567, 0.284231, 0.305546, 0.326471, 0.345985, 0.367151, 
0.386907, 0.408314, 0.429331, 0.448938, 0.470196, 0.490044, 0.509502, 
0.530611, 0.55031, 0.57166, 0.59262, 0.61217, 0.633371, 0.653162, 
0.672563, 0.693616, 0.713258, 0.734551, 0.754434, 0.773927, 0.795071, 
0.814805, 0.83619, 0.857186, 0.876771, 0.898007, 0.917833, 0.93727, 
0.958357, 0.978034, 0.978377, 0.978721, 0.979407, 0.98078, 0.983526, 
0.989017, 0.98936, 0.989704, 0.99039, 0.991763, 0.994509, 0.994852, 
0.995195, 0.995881, 0.997254, 0.997597, 0.997941, 0.998627, 0.99897, 
0.999314, 0.999657, 1.],
y:[2.04082e-8, 0.000306718, 0.000613415, 0.00122681, 0.0024536, 
0.00490718, 0.00981434, 0.0196287, 0.0409084, 0.0607779, 0.0802576, 
0.101388, 0.121109, 0.142481, 0.163463, 0.183035, 0.204257, 0.22407, 
0.243493, 0.264567, 0.284231, 0.305546, 0.326471, 0.345985, 0.367151, 
0.386907, 0.408314, 0.429331, 0.448938, 0.470196, 0.490044, 0.509502, 
0.530611, 0.55031, 0.57166, 0.59262, 0.61217, 0.633371, 0.653162, 
0.672563, 0.693616, 0.713258, 0.734551, 0.754434, 0.773927, 0.795071, 
0.814805, 0.83619, 0.857186, 0.876771, 0.898007, 0.917833, 0.93727, 
0.958357, 0.978034, 0.978377, 0.978721, 0.979407, 0.98078, 0.983526, 
0.989017, 0.98936, 0.989704, 0.99039, 0.991763, 0.994509, 0.994852, 
0.995195, 0.995881, 0.997254, 0.997597, 0.997941, 0.998627, 0.99897,
0.999314, 0.999657, 1.],
mode: 'lines',
name: 'log'
};

var data = [plot0, plot1, plot2, plot3];

var layout = {
  title:  'magnitude of bivector',
  xaxis: { nticks: 10 },
  yaxis: { nticks: 20 },
  height: 400,
  width:  480
};

Plotly.newPlot('fig1', data, layout, {displaylogo: false, autosizable: true});
</script>


<script>

var trace0 = {
x:[2.04082e-8, 0.000306718, 0.000613415, 0.00122681, 0.0024536, 
0.00490718, 0.00981434, 0.0196287, 0.0409084, 0.0607779, 0.0802576, 
0.101388, 0.121109, 0.142481, 0.163463, 0.183035, 0.204257, 0.22407, 
0.243493, 0.264567, 0.284231, 0.305546, 0.326471, 0.345985, 0.367151, 
0.386907, 0.408314, 0.429331, 0.448938, 0.44927, 0.449602, 0.450267, 
0.451595, 0.454253, 0.459567, 0.470196, 0.470506, 0.470816, 0.471437, 
0.472677, 0.475158, 0.48012, 0.48043, 0.48074, 0.48136, 0.482601, 
0.485082, 0.490044, 0.490348, 0.490652, 0.49126, 0.492476, 0.494908, 
0.499773, 0.500077, 0.500381, 0.500989, 0.502205, 0.504637, 0.509502, 
0.509832, 0.510162, 0.510821, 0.512141, 0.514779, 0.515109, 0.515439, 
0.516099, 0.517418, 0.520056, 0.520386, 0.520716, 0.521376, 0.522695, 
0.525334, 0.530611, 0.530919, 0.531227, 0.531842, 0.533073, 0.535536, 
0.535844, 0.536151, 0.536767, 0.537998, 0.54046, 0.540768, 0.541076, 
0.541692, 0.542923, 0.543231, 0.543538, 0.544154, 0.545385, 0.545693, 
0.546001, 0.546616, 0.547848, 0.548155, 0.548463, 0.549079, 0.55031, 
0.550644, 0.550977, 0.551644, 0.551978, 0.552312, 0.552979, 0.553312, 
0.553646, 0.554313, 0.554647, 0.55498, 0.555647, 0.555981, 0.556315, 
0.556648, 0.556982, 0.557315, 0.557649, 0.557983, 0.558316, 0.55865, 
0.558983, 0.559317, 0.559651, 0.559984, 0.560318, 0.560651, 0.560985, 
0.561319, 0.561652, 0.561986, 0.562319, 0.562653, 0.562987, 0.56332, 
0.563654, 0.563987, 0.564321, 0.564654, 0.564988, 0.565322, 0.565655, 
0.566322, 0.566656, 0.56699, 0.567657, 0.56799, 0.568324, 0.568991, 
0.569325, 0.569658, 0.570326, 0.57166, 0.571987, 0.572315, 0.57297, 
0.57428, 0.574607, 0.574935, 0.57559, 0.5769, 0.577228, 0.577555, 
0.57821, 0.57952, 0.58214, 0.582468, 0.582795, 0.58345, 0.58476, 
0.58738, 0.587708, 0.588035, 0.58869, 0.59, 0.59262, 0.592926, 
0.593231, 0.593842, 0.595064, 0.597508, 0.602395, 0.602701, 0.603006, 
0.603617, 0.604839, 0.607283, 0.61217, 0.612501, 0.612833, 0.613495, 
0.61482, 0.617471, 0.622771, 0.623102, 0.623433, 0.624096, 0.625421, 
0.628071, 0.633371, 0.633681, 0.63399, 0.634608, 0.635845, 0.638319, 
0.643267, 0.653162, 0.653465, 0.653769, 0.654375, 0.655587, 0.658013, 
0.662863, 0.672563, 0.672892, 0.673221, 0.673879, 0.675195, 0.677827, 
0.68309, 0.693616, 0.713258, 0.734551, 0.754434, 0.773927, 0.795071, 
0.814805, 0.83619, 0.857186, 0.876771, 0.898007, 0.917833, 0.93727, 
0.958357, 0.978034, 0.978377, 0.978721, 0.979407, 0.98078, 0.983526, 
0.989017, 0.98936, 0.989704, 0.99039, 0.991763, 0.994509, 0.994852, 
0.995195, 0.995881, 0.997254, 0.997597, 0.997941, 0.998627, 0.99897, 
0.999314, 0.999657, 1.],
y:[1.16489e-8, 0.000175073, 0.000350135, 0.000700258, 0.0014005, 
0.00280092, 0.00560138, 0.0111991, 0.0233061, 0.0345468, 0.0454771, 
0.0571998, 0.0679832, 0.0794638, 0.0904918, 0.100531, 0.111113, 
0.120676, 0.129727, 0.139154, 0.147552, 0.156189, 0.164165, 0.171126, 
0.178126, 0.184117, 0.189986, 0.195091, 0.199241, 0.199306, 0.199371, 
0.1995, 0.199756, 0.20026, 0.201235, 0.203044, 0.203094, 0.203144, 
0.203243, 0.20344, 0.203825, 0.204564, 0.204609, 0.204654, 0.204743, 
0.204919, 0.205262, 0.205918, 0.205957, 0.205996, 0.206073, 0.206225, 
0.20652, 0.207082, 0.207115, 0.207149, 0.207215, 0.207347, 0.207601, 
0.20808, 0.208111, 0.208141, 0.208202, 0.208322, 0.208551, 0.208579, 
0.208607, 0.208661, 0.208768, 0.208973, 0.208998, 0.209022, 0.209071, 
0.209165, 0.209345, 0.209666, 0.209683, 0.2097, 0.209733, 0.209798, 
0.20992, 0.209934, 0.209948, 0.209976, 0.21003, 0.210129, 0.21014, 
0.210152, 0.210174, 0.210217, 0.210227, 0.210237, 0.210256, 0.210293, 
0.210302, 0.21031, 0.210327, 0.210358, 0.210366, 0.210373, 0.210387, 
0.210412, 0.210419, 0.210425, 0.210437, 0.210442, 0.210448, 0.210458, 
0.210462, 0.210467, 0.210475, 0.210479, 0.210483, 0.21049, 0.210493, 
0.210496, 0.210498, 0.210501, 0.210503, 0.210505, 0.210507, 0.210508, 
0.21051, 0.210511, 0.210512, 0.210513, 0.210513, 0.210514, 0.210514, 
0.210514, 0.210513, 0.210513, 0.210512, 0.210511, 0.21051, 0.210509, 
0.210507, 0.210505, 0.210503, 0.210501, 0.210498, 0.210496, 0.210493, 
0.21049, 0.210483, 0.210479, 0.210475, 0.210467, 0.210462, 0.210458, 
0.210447, 0.210442, 0.210436, 0.210424, 0.210398, 0.210391, 0.210384, 
0.210369, 0.210336, 0.210328, 0.210319, 0.2103, 0.210261, 0.210251, 
0.21024, 0.210219, 0.210173, 0.210071, 0.210057, 0.210043, 0.210015, 
0.209956, 0.209827, 0.20981, 0.209793, 0.209757, 0.209685, 0.209529, 
0.20951, 0.209491, 0.209452, 0.209371, 0.209202, 0.208828, 0.208803, 
0.208777, 0.208726, 0.208622, 0.208405, 0.207935, 0.207901, 0.207867, 
0.207799, 0.20766, 0.20737, 0.206748, 0.206708, 0.206667, 0.206584, 
0.206416, 0.206069, 0.205332, 0.205287, 0.205242, 0.205151, 0.204968, 
0.204591, 0.2038, 0.202063, 0.202006, 0.20195, 0.201836, 0.201606, 
0.201136, 0.20016, 0.198056, 0.197981, 0.197906, 0.197755, 0.197451, 
0.19683, 0.195545, 0.192793, 0.187009, 0.179771, 0.172089, 0.16368, 
0.153564, 0.14318, 0.130887, 0.117757, 0.104553, 0.0891866, 0.073849, 
0.0578795, 0.0395043, 0.0213705, 0.0210458, 0.0207207, 0.0200697, 
0.0187643, 0.0161395, 0.0108341, 0.0105, 0.0101657, 0.00949609, 
0.00815346, 0.00545424, 0.00511554, 0.00477654, 0.00409766, 
0.00273643, 0.0023954, 0.00205407, 0.00137055, 0.00102836, 
0.000685868, 0.000343089, 2.04082e-8],
mode: 'lines',
name: 'standard'
};

var trace1 = {
x:[2.04082e-8, 0.000306718, 0.000613415, 0.00122681, 0.0024536, 
0.00490718, 0.00981434, 0.0196287, 0.0409084, 0.0607779, 0.0802576, 
0.101388, 0.121109, 0.142481, 0.163463, 0.183035, 0.204257, 0.22407, 
0.243493, 0.264567, 0.284231, 0.305546, 0.326471, 0.345985, 0.367151, 
0.386907, 0.408314, 0.429331, 0.448938, 0.44927, 0.449602, 0.450267, 
0.451595, 0.454253, 0.459567, 0.470196, 0.470506, 0.470816, 0.471437, 
0.472677, 0.475158, 0.48012, 0.490044, 0.490348, 0.490652, 0.49126, 
0.492476, 0.494908, 0.499773, 0.500077, 0.500381, 0.500989, 0.502205, 
0.504637, 0.509502, 0.509832, 0.510162, 0.510821, 0.512141, 0.514779, 
0.520056, 0.520386, 0.520716, 0.521376, 0.522695, 0.525334, 0.530611, 
0.530919, 0.531227, 0.531842, 0.533073, 0.535536, 0.535844, 0.536151, 
0.536767, 0.537998, 0.54046, 0.540768, 0.541076, 0.541692, 0.542923, 
0.545385, 0.545693, 0.546001, 0.546616, 0.547848, 0.55031, 0.550644, 
0.550977, 0.551644, 0.552979, 0.553312, 0.553646, 0.554313, 0.555647, 
0.555981, 0.556315, 0.556982, 0.558316, 0.55865, 0.558983, 0.559651, 
0.560985, 0.561319, 0.561652, 0.562319, 0.562653, 0.562987, 0.563654, 
0.563987, 0.564321, 0.564988, 0.565322, 0.565655, 0.566322, 0.566656, 
0.56699, 0.567323, 0.567657, 0.56799, 0.568324, 0.568991, 0.569325, 
0.569658, 0.569992, 0.570326, 0.570659, 0.570993, 0.571326, 0.57166, 
0.571987, 0.572315, 0.572642, 0.57297, 0.573297, 0.573625, 0.573952, 
0.57428, 0.574607, 0.574935, 0.575262, 0.57559, 0.575917, 0.576245, 
0.576572, 0.5769, 0.577228, 0.577555, 0.577883, 0.57821, 0.578538, 
0.578865, 0.57952, 0.579848, 0.580175, 0.58083, 0.581158, 0.581485, 
0.58214, 0.582468, 0.582795, 0.58345, 0.58476, 0.585088, 0.585415, 
0.58607, 0.58738, 0.587708, 0.588035, 0.58869, 0.59, 0.59262, 
0.592926, 0.593231, 0.593842, 0.595064, 0.597508, 0.597813, 0.598119, 
0.59873, 0.599951, 0.602395, 0.602701, 0.603006, 0.603617, 0.604839, 
0.607283, 0.61217, 0.612501, 0.612833, 0.613495, 0.61482, 0.617471, 
0.622771, 0.623102, 0.623433, 0.624096, 0.625421, 0.628071, 0.633371, 
0.633681, 0.63399, 0.634608, 0.635845, 0.638319, 0.643267, 0.643576, 
0.643885, 0.644504, 0.645741, 0.648215, 0.653162, 0.653465, 0.653769, 
0.654375, 0.655587, 0.658013, 0.662863, 0.672563, 0.672892, 0.673221, 
0.673879, 0.675195, 0.677827, 0.68309, 0.693616, 0.693923, 0.694229, 
0.694843, 0.696071, 0.698526, 0.703437, 0.713258, 0.734551, 0.754434, 
0.773927, 0.795071, 0.814805, 0.83619, 0.857186, 0.876771, 0.898007, 
0.917833, 0.93727, 0.958357, 0.978034, 0.978377, 0.978721, 0.979407, 
0.98078, 0.983526, 0.989017, 0.98936, 0.989704, 0.99039, 0.991763, 
0.994509, 0.994852, 0.995195, 0.995881, 0.997254, 0.997597, 0.997941, 
0.998627, 0.99897, 0.999314, 0.999657, 1.],
y:[2.25961e-9, 0.00003396, 0.0000679178, 0.000135833, 0.000271663, 
0.000543313, 0.00108654, 0.00217244, 0.00452159, 0.00670374, 
0.00882716, 0.0111068, 0.0132065, 0.0154455, 0.0176004, 0.0195662, 
0.0216437, 0.0235265, 0.0253142, 0.0271829, 0.0288547, 0.0305823, 
0.0321867, 0.0335958, 0.0350232, 0.0362553, 0.0374752, 0.0385504, 
0.0394386, 0.0394527, 0.0394667, 0.0394946, 0.0395502, 0.0396596, 
0.039872, 0.0402706, 0.0402817, 0.0402928, 0.0403148, 0.0403585, 
0.0404445, 0.0406107, 0.0409191, 0.040928, 0.0409369, 0.0409547, 
0.0409898, 0.0410585, 0.0411901, 0.0411981, 0.041206, 0.0412218, 
0.041253, 0.0413139, 0.0414296, 0.0414372, 0.0414447, 0.0414596, 
0.0414891, 0.0415461, 0.041653, 0.0416594, 0.0416657, 0.0416783, 
0.0417029, 0.0417503, 0.0418378, 0.0418426, 0.0418474, 0.0418568, 
0.0418753, 0.0419106, 0.0419149, 0.0419191, 0.0419275, 0.0419438, 
0.0419748, 0.0419785, 0.0419822, 0.0419895, 0.0420036, 0.0420303, 
0.0420335, 0.0420366, 0.0420428, 0.0420547, 0.042077, 0.0420798, 
0.0420826, 0.0420881, 0.0420986, 0.0421011, 0.0421035, 0.0421084, 
0.0421175, 0.0421197, 0.0421219, 0.042126, 0.0421339, 0.0421357, 
0.0421376, 0.0421411, 0.0421476, 0.0421491, 0.0421506, 0.0421535, 
0.0421548, 0.0421562, 0.0421587, 0.0421599, 0.042161, 0.0421632, 
0.0421642, 0.0421652, 0.0421671, 0.0421679, 0.0421688, 0.0421695, 
0.0421703, 0.042171, 0.0421716, 0.0421728, 0.0421733, 0.0421738, 
0.0421743, 0.0421747, 0.042175, 0.0421753, 0.0421756, 0.0421759, 
0.042176, 0.0421762, 0.0421763, 0.0421764, 0.0421764, 0.0421764, 
0.0421763, 0.0421762, 0.0421761, 0.0421759, 0.0421757, 0.0421754, 
0.0421751, 0.0421748, 0.0421744, 0.042174, 0.0421735, 0.042173, 
0.0421725, 0.0421719, 0.0421712, 0.0421706, 0.0421691, 0.0421683, 
0.0421675, 0.0421657, 0.0421647, 0.0421637, 0.0421616, 0.0421605, 
0.0421593, 0.0421569, 0.0421515, 0.04215, 0.0421485, 0.0421454, 
0.0421387, 0.0421369, 0.0421351, 0.0421313, 0.0421232, 0.0421051, 
0.0421028, 0.0421005, 0.0420957, 0.0420858, 0.0420641, 0.0420612, 
0.0420583, 0.0420523, 0.04204, 0.0420136, 0.0420101, 0.0420066, 
0.0419995, 0.0419848, 0.0419537, 0.0418842, 0.0418791, 0.041874, 
0.0418637, 0.0418424, 0.0417979, 0.0417002, 0.0416937, 0.0416872, 
0.041674, 0.041647, 0.041591, 0.0414701, 0.0414627, 0.0414553, 
0.0414403, 0.0414098, 0.0413468, 0.0412132, 0.0412045, 0.0411958, 
0.0411782, 0.0411425, 0.0410693, 0.0409149, 0.0409051, 0.0408953, 
0.0408755, 0.0408355, 0.0407535, 0.0405818, 0.0402079, 0.0401945, 
0.0401811, 0.040154, 0.0400994, 0.0399878, 0.0397553, 0.0392534, 
0.039238, 0.0392226, 0.0391916, 0.0391291, 0.0390021, 0.0387399, 
0.0381822, 0.0368189, 0.0353511, 0.0337252, 0.0317475, 0.029696, 
0.0272437, 0.0245992, 0.0219163, 0.0187668, 0.0155963, 0.0122685, 
0.00840768, 0.00456593, 0.00449685, 0.0044277, 0.00428918, 
0.00401129, 0.00345207, 0.00231985, 0.00224847, 0.00217703, 
0.00203392, 0.00174683, 0.00116919, 0.00109666, 0.00102405, 
0.000878631, 0.000586916, 0.000513806, 0.000440624, 0.000294041, 
0.000220641, 0.000147168, 0.0000736225, 4.37963e-9],
mode: 'lines',
name: 'half angle'
};

var trace2 = {
x:[2.04082e-8, 0.000306718, 0.000613415, 0.00122681, 0.0024536, 
0.00490718, 0.00981434, 0.0196287, 0.0409084, 0.0607779, 0.0802576, 
0.101388, 0.121109, 0.142481, 0.163463, 0.183035, 0.204257, 0.22407, 
0.243493, 0.264567, 0.284231, 0.305546, 0.326471, 0.345985, 0.367151, 
0.386907, 0.408314, 0.429331, 0.448938, 0.470196, 0.470506, 0.470816, 
0.471437, 0.472677, 0.475158, 0.48012, 0.490044, 0.490348, 0.490652, 
0.49126, 0.492476, 0.494908, 0.499773, 0.500077, 0.500381, 0.500989, 
0.502205, 0.504637, 0.509502, 0.509832, 0.510162, 0.510821, 0.512141, 
0.514779, 0.520056, 0.520386, 0.520716, 0.521376, 0.522695, 0.525334, 
0.530611, 0.530919, 0.531227, 0.531842, 0.533073, 0.535536, 0.535844, 
0.536151, 0.536767, 0.537998, 0.54046, 0.540768, 0.541076, 0.541692, 
0.542923, 0.545385, 0.545693, 0.546001, 0.546616, 0.547848, 0.55031, 
0.550644, 0.550977, 0.551644, 0.552979, 0.553312, 0.553646, 0.554313, 
0.555647, 0.555981, 0.556315, 0.556982, 0.558316, 0.560985, 0.561319, 
0.561652, 0.562319, 0.563654, 0.563987, 0.564321, 0.564988, 0.566322, 
0.566656, 0.56699, 0.567657, 0.56799, 0.568324, 0.568991, 0.569325, 
0.569658, 0.570326, 0.570659, 0.570993, 0.57166, 0.571987, 0.572315, 
0.572642, 0.57297, 0.573297, 0.573625, 0.573952, 0.57428, 0.574607, 
0.574935, 0.575262, 0.57559, 0.575917, 0.576245, 0.576572, 0.5769, 
0.577228, 0.577555, 0.577883, 0.57821, 0.578538, 0.578865, 0.579193, 
0.57952, 0.579848, 0.580175, 0.580503, 0.58083, 0.581158, 0.581485, 
0.58214, 0.582468, 0.582795, 0.58345, 0.583778, 0.584105, 0.58476, 
0.585088, 0.585415, 0.58607, 0.58738, 0.587708, 0.588035, 0.58869, 
0.59, 0.590328, 0.590655, 0.59131, 0.59262, 0.592926, 0.593231, 
0.593842, 0.595064, 0.597508, 0.597813, 0.598119, 0.59873, 0.599951, 
0.602395, 0.602701, 0.603006, 0.603617, 0.604839, 0.607283, 0.61217, 
0.612501, 0.612833, 0.613495, 0.61482, 0.617471, 0.622771, 0.623102, 
0.623433, 0.624096, 0.625421, 0.628071, 0.633371, 0.633681, 0.63399, 
0.634608, 0.635845, 0.638319, 0.643267, 0.643576, 0.643885, 0.644504, 
0.645741, 0.648215, 0.653162, 0.653465, 0.653769, 0.654375, 0.655587, 
0.658013, 0.662863, 0.672563, 0.672892, 0.673221, 0.673879, 0.675195, 
0.677827, 0.68309, 0.693616, 0.693923, 0.694229, 0.694843, 0.696071, 
0.698526, 0.703437, 0.713258, 0.734551, 0.754434, 0.773927, 0.795071, 
0.814805, 0.83619, 0.857186, 0.876771, 0.898007, 0.917833, 0.93727, 
0.958357, 0.978034, 0.978377, 0.978721, 0.979407, 0.98078, 0.983526, 
0.989017, 0.98936, 0.989704, 0.99039, 0.991763, 0.994509, 0.994852, 
0.995195, 0.995881, 0.997254, 0.997597, 0.997941, 0.998627, 0.99897, 
0.999314, 0.999657, 1.],
y:[5.34126e-10, 8.02747e-6, 0.0000160544, 0.0000321082, 0.0000642156, 
0.000128428, 0.000256837, 0.000513525, 0.00106885, 0.00158477, 
0.00208688, 0.00262607, 0.00312284, 0.00365275, 0.004163, 0.00462872, 
0.00512117, 0.00556779, 0.00599215, 0.00643613, 0.00683368, 
0.00724498, 0.00762745, 0.00796384, 0.00830516, 0.00860036, 
0.00889332, 0.00915228, 0.00936696, 0.00956898, 0.00957168, 
0.00957438, 0.00957975, 0.0095904, 0.00961137, 0.00965193, 
0.00972745, 0.00972964, 0.00973183, 0.00973618, 0.0097448, 0.0097617, 
0.00979411, 0.00979608, 0.00979803, 0.00980192, 0.00980962, 
0.00982466, 0.00985334, 0.00985521, 0.00985708, 0.00986079, 
0.0098681, 0.0098823, 0.00990901, 0.00991061, 0.00991219, 0.00991534, 
0.00992152, 0.00993345, 0.00995558, 0.0099568, 0.00995802, 
0.00996042, 0.00996512, 0.00997415, 0.00997524, 0.00997633, 
0.00997847, 0.00998267, 0.00999067, 0.00999164, 0.00999259, 
0.00999448, 0.00999816, 0.0100051, 0.010006, 0.0100068, 0.0100084, 
0.0100116, 0.0100175, 0.0100183, 0.010019, 0.0100205, 0.0100234, 
0.010024, 0.0100247, 0.010026, 0.0100286, 0.0100292, 0.0100298, 
0.010031, 0.0100332, 0.0100372, 0.0100376, 0.0100381, 0.0100389, 
0.0100405, 0.0100409, 0.0100412, 0.0100419, 0.0100432, 0.0100435, 
0.0100438, 0.0100443, 0.0100446, 0.0100448, 0.0100453, 0.0100455, 
0.0100457, 0.0100461, 0.0100463, 0.0100464, 0.0100468, 0.0100469, 
0.010047, 0.0100471, 0.0100472, 0.0100473, 0.0100474, 0.0100475, 
0.0100476, 0.0100476, 0.0100477, 0.0100477, 0.0100477, 0.0100477, 
0.0100477, 0.0100477, 0.0100477, 0.0100477, 0.0100477, 0.0100476, 
0.0100476, 0.0100475, 0.0100475, 0.0100474, 0.0100473, 0.0100472, 
0.0100471, 0.010047, 0.0100468, 0.0100467, 0.0100466, 0.0100462, 
0.0100461, 0.0100459, 0.0100455, 0.0100452, 0.010045, 0.0100445, 
0.0100443, 0.010044, 0.0100435, 0.0100422, 0.0100419, 0.0100415, 
0.0100408, 0.0100393, 0.0100389, 0.0100384, 0.0100376, 0.0100357, 
0.0100352, 0.0100347, 0.0100338, 0.0100318, 0.0100273, 0.0100267, 
0.0100261, 0.0100248, 0.0100222, 0.0100166, 0.0100159, 0.0100152, 
0.0100136, 0.0100105, 0.0100037, 0.00998858, 0.00998747, 0.00998634, 
0.00998407, 0.00997939, 0.00996954, 0.00994778, 0.00994633, 
0.00994487, 0.00994191, 0.00993587, 0.00992328, 0.00989601, 
0.00989434, 0.00989265, 0.00988925, 0.00988234, 0.00986805, 
0.00983764, 0.00983566, 0.00983366, 0.00982965, 0.00982151, 
0.00980476, 0.0097694, 0.00976715, 0.00976489, 0.00976035, 
0.00975115, 0.0097323, 0.00969279, 0.00960642, 0.00960332, 0.0096002, 
0.00959394, 0.00958129, 0.00955542, 0.00950149, 0.00938474, 
0.00938115, 0.00937756, 0.00937034, 0.00935578, 0.00932617, 
0.00926497, 0.0091346, 0.00881486, 0.00846939, 0.00808555, 
0.00761736, 0.00713049, 0.00654712, 0.00591657, 0.00527551, 
0.00452142, 0.00376077, 0.00296086, 0.00203101, 0.00110397, 
0.00108728, 0.00107058, 0.00103712, 0.000969988, 0.000834867, 
0.000561188, 0.000543931, 0.000526655, 0.00049205, 0.000422625, 
0.000282908, 0.000265362, 0.000247798, 0.000212616, 0.000142034, 
0.000124343, 0.000106635, 0.0000711627, 0.0000533995, 0.0000356182, 
0.0000178187, 1.06001e-9],
mode: 'lines',
name: 'half angle/half angle'

};

var trace3 = {
x:[2.04082e-8, 0.000306718, 0.000613415, 0.00122681, 0.0024536, 
0.00490718, 0.00981434, 0.0196287, 0.0409084, 0.0607779, 0.0802576, 
0.101388, 0.121109, 0.142481, 0.163463, 0.183035, 0.204257, 0.22407, 
0.243493, 0.264567, 0.284231, 0.305546, 0.326471, 0.345985, 0.367151, 
0.386907, 0.408314, 0.429331, 0.448938, 0.470196, 0.490044, 0.509502, 
0.509832, 0.510162, 0.510821, 0.512141, 0.514779, 0.520056, 0.530611, 
0.530919, 0.531227, 0.531842, 0.533073, 0.535536, 0.54046, 0.540768, 
0.541076, 0.541692, 0.542923, 0.545385, 0.55031, 0.550644, 0.550977, 
0.551644, 0.552979, 0.555647, 0.560985, 0.561319, 0.561652, 0.562319, 
0.563654, 0.566322, 0.57166, 0.571987, 0.572315, 0.57297, 0.57428, 
0.5769, 0.577228, 0.577555, 0.57821, 0.57952, 0.58214, 0.582468, 
0.582795, 0.58345, 0.58476, 0.58738, 0.587708, 0.588035, 0.58869, 
0.59, 0.59262, 0.592926, 0.593231, 0.593842, 0.595064, 0.595369, 
0.595675, 0.596286, 0.597508, 0.597813, 0.598119, 0.59873, 0.599951, 
0.600257, 0.600562, 0.601173, 0.602395, 0.602701, 0.603006, 0.603617, 
0.603923, 0.604228, 0.604839, 0.605144, 0.60545, 0.606061, 0.606366, 
0.606672, 0.607283, 0.607588, 0.607894, 0.608199, 0.608505, 0.60881, 
0.609116, 0.609726, 0.610032, 0.610337, 0.610643, 0.610948, 0.611254, 
0.611559, 0.611865, 0.61217, 0.612501, 0.612833, 0.613164, 0.613495, 
0.613827, 0.614158, 0.614489, 0.61482, 0.615152, 0.615483, 0.615814, 
0.616145, 0.616477, 0.616808, 0.617139, 0.617471, 0.617802, 0.618133, 
0.618796, 0.619127, 0.619458, 0.620121, 0.620452, 0.620783, 0.621446, 
0.622771, 0.623102, 0.623433, 0.624096, 0.625421, 0.625752, 0.626083, 
0.626746, 0.628071, 0.628402, 0.628734, 0.629396, 0.630721, 0.633371, 
0.633681, 0.63399, 0.634608, 0.635845, 0.638319, 0.638628, 0.638938, 
0.639556, 0.640793, 0.643267, 0.643576, 0.643885, 0.644504, 0.645741, 
0.648215, 0.653162, 0.653465, 0.653769, 0.654375, 0.655587, 0.658013, 
0.662863, 0.663166, 0.663469, 0.664075, 0.665288, 0.667713, 0.672563, 
0.672892, 0.673221, 0.673879, 0.675195, 0.677827, 0.68309, 0.683419, 
0.683747, 0.684405, 0.685721, 0.688353, 0.693616, 0.693923, 0.694229, 
0.694843, 0.696071, 0.698526, 0.703437, 0.713258, 0.71359, 0.713923, 
0.714589, 0.715919, 0.718581, 0.723904, 0.734551, 0.734862, 0.735172, 
0.735794, 0.737036, 0.739522, 0.744492, 0.754434, 0.773927, 0.795071, 
0.814805, 0.83619, 0.857186, 0.876771, 0.898007, 0.917833, 0.93727, 
0.958357, 0.978034, 0.978377, 0.978721, 0.979407, 0.98078, 0.983526, 
0.989017, 0.98936, 0.989704, 0.99039, 0.991763, 0.994509, 0.994852, 
0.995195, 0.995881, 0.997254, 0.997597, 0.997941, 0.998627, 0.99897, 
0.999314, 0.999657, 1.],
y:[-4.37963e-9, -0.0000658222, -0.00013164, -0.000263275, 
-0.000526545, -0.00105307, -0.00210602, -0.00421113, -0.00876795, 
-0.0130068, -0.0171398, -0.0215894, -0.0257023, -0.0301072, 
-0.0343694, -0.038281, -0.0424435, -0.0462463, -0.0498881, 
-0.0537334, -0.0572128, -0.0608554, -0.0642901, -0.0673571, 
-0.0705238, -0.0733187, -0.0761597, -0.0787461, -0.0809659, 
-0.0831482, -0.0849621, -0.0865179, -0.0865423, -0.0865666, 
-0.0866151, -0.0867112, -0.0869002, -0.0872652, -0.0879418, 
-0.0879605, -0.0879791, -0.0880161, -0.0880893, -0.0882328, 
-0.0885077, -0.0885243, -0.0885409, -0.0885739, -0.0886391, 
-0.0887663, -0.0890084, -0.0890242, -0.08904, -0.0890712, -0.0891327, 
-0.0892519, -0.0894755, -0.0894888, -0.089502, -0.0895282, 
-0.0895796, -0.0896787, -0.0898614, -0.089872, -0.0898824, 
-0.0899031, -0.0899435, -0.0900205, -0.0900298, -0.090039, 
-0.0900571, -0.0900924, -0.0901592, -0.0901672, -0.0901751, 
-0.0901907, -0.0902208, -0.0902772, -0.0902839, -0.0902905, 
-0.0903034, -0.0903283, -0.0903742, -0.0903792, -0.0903841, 
-0.0903937, -0.0904121, -0.0904166, -0.0904209, -0.0904294, 
-0.0904454, -0.0904493, -0.090453, -0.0904603, -0.090474, -0.0904773, 
-0.0904804, -0.0904866, -0.0904979, -0.0905006, -0.0905031, 
-0.0905081, -0.0905104, -0.0905127, -0.090517, -0.0905191, -0.090521, 
-0.0905248, -0.0905265, -0.0905282, -0.0905313, -0.0905328, 
-0.0905342, -0.0905355, -0.0905367, -0.0905378, -0.0905389, 
-0.0905408, -0.0905417, -0.0905424, -0.0905431, -0.0905438, 
-0.0905443, -0.0905448, -0.0905452, -0.0905455, -0.0905457, 
-0.0905459, -0.090546, -0.0905459, -0.0905458, -0.0905456, 
-0.0905454, -0.090545, -0.0905445, -0.090544, -0.0905433, -0.0905426, 
-0.0905417, -0.0905408, -0.0905398, -0.0905387, -0.0905375, 
-0.0905362, -0.0905334, -0.0905318, -0.0905302, -0.0905266, 
-0.0905246, -0.0905226, -0.0905183, -0.0905086, -0.0905059, 
-0.0905032, -0.0904974, -0.0904847, -0.0904813, -0.0904778, 
-0.0904705, -0.0904548, -0.0904506, -0.0904464, -0.0904376, 
-0.0904189, -0.090377, -0.0903717, -0.0903663, -0.0903553, 
-0.0903324, -0.0902824, -0.0902758, -0.090269, -0.0902554, -0.090227, 
-0.0901662, -0.0901583, -0.0901502, -0.0901338, -0.0901, -0.0900283, 
-0.0898682, -0.0898577, -0.0898471, -0.0898256, -0.0897816, 
-0.0896896, -0.0894891, -0.0894758, -0.0894625, -0.0894355, 
-0.0893806, -0.0892665, -0.0890215, -0.089004, -0.0889865, 
-0.0889511, -0.0888791, -0.0887299, -0.0884113, -0.0883904, 
-0.0883695, -0.0883273, -0.0882417, -0.0880651, -0.087691, 
-0.0876683, -0.0876456, -0.0875997, -0.0875069, -0.0873165, 
-0.086917, -0.0860417, -0.0860102, -0.0859786, -0.0859151, 
-0.0857866, -0.0855239, -0.0849752, -0.0837829, -0.0837462, 
-0.0837094, -0.0836354, -0.0834861, -0.0831821, -0.0825527, 
-0.0812064, -0.0782172, -0.0744266, -0.0703445, -0.0652921, 
-0.0596572, -0.0537618, -0.0466297, -0.0392345, -0.0312519, 
-0.0217199, -0.0119561, -0.011778, -0.0115997, -0.0112422, -0.010524, 
-0.00907433, -0.00612187, -0.00593496, -0.00574778, -0.00537256, 
-0.00461875, -0.00309751, -0.00290607, -0.00271435, -0.00233004, 
-0.00155798, -0.00136425, -0.00117022, -0.000781311, -0.000586421, 
-0.000391242, -0.000195772, -1.16489e-8],
mode: 'lines',
name: 'Cayley'

};

var trace4 = {
x:[2.04082e-8, 0.000306718, 0.000613415, 0.00122681, 0.0024536, 
0.00490718, 0.00981434, 0.0196287, 0.0409084, 0.0607779, 0.0802576, 
0.101388, 0.121109, 0.142481, 0.163463, 0.183035, 0.204257, 0.22407, 
0.243493, 0.264567, 0.284231, 0.305546, 0.326471, 0.345985, 0.367151, 
0.386907, 0.408314, 0.429331, 0.448938, 0.470196, 0.470506, 0.470816, 
0.471437, 0.472677, 0.475158, 0.48012, 0.490044, 0.490348, 0.490652, 
0.49126, 0.492476, 0.494908, 0.499773, 0.500077, 0.500381, 0.500989, 
0.502205, 0.504637, 0.509502, 0.509832, 0.510162, 0.510821, 0.512141, 
0.514779, 0.520056, 0.520386, 0.520716, 0.521376, 0.522695, 0.525334, 
0.530611, 0.530919, 0.531227, 0.531842, 0.533073, 0.535536, 0.54046, 
0.540768, 0.541076, 0.541692, 0.542923, 0.545385, 0.55031, 0.550644, 
0.550977, 0.551644, 0.552979, 0.555647, 0.555981, 0.556315, 0.556982, 
0.558316, 0.560985, 0.561319, 0.561652, 0.562319, 0.563654, 0.563987, 
0.564321, 0.564988, 0.566322, 0.566656, 0.56699, 0.567657, 0.568991, 
0.569325, 0.569658, 0.570326, 0.57166, 0.571987, 0.572315, 0.57297, 
0.573297, 0.573625, 0.57428, 0.574607, 0.574935, 0.57559, 0.5769, 
0.577228, 0.577555, 0.57821, 0.578538, 0.578865, 0.57952, 0.579848, 
0.580175, 0.580503, 0.58083, 0.581158, 0.581485, 0.581813, 0.58214, 
0.582468, 0.582795, 0.583123, 0.58345, 0.583778, 0.584105, 0.584433, 
0.58476, 0.585088, 0.585415, 0.585743, 0.58607, 0.586398, 0.586725, 
0.587053, 0.58738, 0.587708, 0.588035, 0.588363, 0.58869, 0.589018, 
0.589345, 0.59, 0.590328, 0.590655, 0.59131, 0.591638, 0.591965, 
0.59262, 0.592926, 0.593231, 0.593842, 0.594148, 0.594453, 0.595064, 
0.595369, 0.595675, 0.596286, 0.597508, 0.597813, 0.598119, 0.59873, 
0.599951, 0.602395, 0.602701, 0.603006, 0.603617, 0.604839, 0.607283, 
0.607588, 0.607894, 0.608505, 0.609726, 0.61217, 0.612501, 0.612833, 
0.613495, 0.61482, 0.617471, 0.617802, 0.618133, 0.618796, 0.620121, 
0.622771, 0.623102, 0.623433, 0.624096, 0.625421, 0.628071, 0.633371, 
0.633681, 0.63399, 0.634608, 0.635845, 0.638319, 0.643267, 0.643576, 
0.643885, 0.644504, 0.645741, 0.648215, 0.653162, 0.653465, 0.653769, 
0.654375, 0.655587, 0.658013, 0.662863, 0.663166, 0.663469, 0.664075, 
0.665288, 0.667713, 0.672563, 0.672892, 0.673221, 0.673879, 0.675195, 
0.677827, 0.68309, 0.693616, 0.693923, 0.694229, 0.694843, 0.696071, 
0.698526, 0.703437, 0.713258, 0.71359, 0.713923, 0.714589, 0.715919, 
0.718581, 0.723904, 0.734551, 0.754434, 0.773927, 0.795071, 0.814805, 
0.83619, 0.857186, 0.876771, 0.898007, 0.917833, 0.93727, 0.958357, 
0.978034, 0.978377, 0.978721, 0.979407, 0.98078, 0.983526, 0.989017, 
0.98936, 0.989704, 0.99039, 0.991763, 0.994509, 0.994852, 0.995195, 
0.995881, 0.997254, 0.997597, 0.997941, 0.998627, 0.99897, 0.999314, 
0.999657, 1.],
y:[-1.06001e-9, -0.0000159311, -0.0000318611, -0.0000637211, 
-0.000127441, -0.000254876, -0.000509716, -0.00101915, -0.00212147, 
-0.00314589, -0.00414342, -0.00521535, -0.00620383, -0.00725939, 
-0.00827714, -0.00920746, -0.0101929, -0.0110884, -0.011941, 
-0.0128353, -0.0136384, -0.014472, -0.0152501, -0.0159373, 
-0.0166379, -0.0172472, -0.017856, -0.0183987, -0.018853, -0.0192861, 
-0.019292, -0.0192978, -0.0193094, -0.0193325, -0.019378, -0.0194662, 
-0.0196318, -0.0196367, -0.0196415, -0.0196511, -0.0196702, 
-0.0197076, -0.0197798, -0.0197842, -0.0197886, -0.0197973, 
-0.0198145, -0.0198483, -0.0199132, -0.0199174, -0.0199217, 
-0.0199301, -0.0199468, -0.0199793, -0.020041, -0.0200448, 
-0.0200484, -0.0200558, -0.0200702, -0.0200983, -0.020151, 
-0.0201539, -0.0201568, -0.0201626, -0.020174, -0.0201961, 
-0.0202371, -0.0202395, -0.0202419, -0.0202467, -0.0202561, 
-0.020274, -0.0203069, -0.020309, -0.020311, -0.0203151, -0.0203229, 
-0.0203378, -0.0203396, -0.0203413, -0.0203447, -0.0203514, 
-0.0203638, -0.0203652, -0.0203666, -0.0203695, -0.0203749, 
-0.0203762, -0.0203775, -0.02038, -0.0203847, -0.0203859, -0.020387, 
-0.0203892, -0.0203933, -0.0203943, -0.0203953, -0.0203972, 
-0.0204007, -0.0204015, -0.0204023, -0.0204038, -0.0204045, 
-0.0204053, -0.0204066, -0.0204073, -0.0204079, -0.0204091, 
-0.0204114, -0.0204119, -0.0204123, -0.0204133, -0.0204137, 
-0.0204141, -0.0204148, -0.0204152, -0.0204155, -0.0204158, 
-0.0204161, -0.0204164, -0.0204166, -0.0204169, -0.0204171, 
-0.0204173, -0.0204174, -0.0204176, -0.0204177, -0.0204178, 
-0.0204179, -0.020418, -0.0204181, -0.0204181, -0.0204181, 
-0.0204181, -0.0204181, -0.020418, -0.020418, -0.0204179, -0.0204178, 
-0.0204177, -0.0204175, -0.0204174, -0.0204172, -0.020417, 
-0.0204168, -0.0204162, -0.020416, -0.0204157, -0.020415, -0.0204146, 
-0.0204143, -0.0204134, -0.020413, -0.0204126, -0.0204117, 
-0.0204112, -0.0204107, -0.0204097, -0.0204091, -0.0204085, 
-0.0204073, -0.0204048, -0.0204041, -0.0204034, -0.0204019, 
-0.0203987, -0.0203916, -0.0203906, -0.0203896, -0.0203876, 
-0.0203833, -0.0203739, -0.0203726, -0.0203713, -0.0203687, 
-0.0203633, -0.0203516, -0.0203499, -0.0203482, -0.0203448, 
-0.0203376, -0.0203222, -0.0203202, -0.0203181, -0.020314, 
-0.0203054, -0.0202873, -0.0202849, -0.0202826, -0.0202777, 
-0.0202678, -0.0202469, -0.0202008, -0.020198, -0.0201951, 
-0.0201893, -0.0201774, -0.0201528, -0.0200997, -0.0200962, 
-0.0200927, -0.0200857, -0.0200713, -0.0200417, -0.0199785, 
-0.0199745, -0.0199704, -0.0199623, -0.0199457, -0.0199117, 
-0.0198399, -0.0198352, -0.0198306, -0.0198212, -0.0198021, 
-0.0197631, -0.0196813, -0.0196755, -0.0196698, -0.0196582, 
-0.0196347, -0.0195867, -0.0194861, -0.0192666, -0.0192598, 
-0.019253, -0.0192394, -0.0192118, -0.0191557, -0.0190394, 
-0.0187901, -0.0187813, -0.0187725, -0.0187547, -0.0187188, 
-0.0186457, -0.0184947, -0.0181725, -0.0174974, -0.0167402, 
-0.0158085, -0.0148321, -0.0136534, -0.0123702, -0.0110571, 
-0.0095026, -0.00792472, -0.00625558, -0.0043036, -0.0023458, 
-0.00231046, -0.00227508, -0.00220419, -0.00206192, -0.0017754, 
-0.00119435, -0.00115768, -0.00112097, -0.00104742, -0.000899814, 
-0.000602583, -0.000565239, -0.000527854, -0.000452955, -0.000302649, 
-0.000264967, -0.000227242, -0.000151666, -0.000113814, -0.000075919,
-0.0000379819, -2.25961e-9],
mode: 'lines',
name: 'half angle/Cayley'

};

var trace5 = { x:[2.04082e-8, 1.], y:[0., 0.], mode: 'lines', name: 'log' };


var data = [trace0, trace1, trace3, trace5, trace2, trace4];

var layout = {
  title:  'magnitude of bivector vs. linear',
  xaxis: { nticks: 10 },
  yaxis: { nticks: 12 },
  height: 400,
  width:  600
};

Plotly.newPlot('fig2', data, layout, {displaylogo: false, autosizable: true});
</script>


<script>
var plot0 = {
x:[2.04082e-8, 0.000306718, 0.000613415, 0.000920113, 0.00122681, 
0.00153351, 0.00184021, 0.0021469, 0.0024536, 0.0027603, 0.003067, 
0.00337369, 0.00368039, 0.00398709, 0.00429379, 0.00490718, 
0.00521388, 0.00552058, 0.00613397, 0.00644067, 0.00674737, 
0.00736076, 0.00766746, 0.00797416, 0.00858755, 0.00981434, 0.010121, 
0.0104277, 0.0110411, 0.0122679, 0.0125746, 0.0128813, 0.0134947, 
0.0147215, 0.0150282, 0.0153349, 0.0159483, 0.0171751, 0.0196287, 
0.0199612, 0.0202937, 0.0209586, 0.0222886, 0.0249486, 0.0252811, 
0.0256136, 0.0262786, 0.0276085, 0.0302685, 0.030601, 0.0309335, 
0.0315985, 0.0329285, 0.0355884, 0.0409084, 0.0412188, 0.0415293, 
0.0421502, 0.043392, 0.0458757, 0.0508431, 0.0511536, 0.051464, 
0.052085, 0.0533268, 0.0558105, 0.0607779, 0.0610823, 0.0613866, 
0.0619954, 0.0632129, 0.0656478, 0.0705178, 0.0708221, 0.0711265, 
0.0717352, 0.0729527, 0.0753877, 0.0802576, 0.0805878, 0.080918, 
0.0815783, 0.082899, 0.0855404, 0.0908231, 0.101388, 0.101697, 
0.102005, 0.102621, 0.103854, 0.106319, 0.111249, 0.121109, 0.142481, 
0.163463, 0.183035, 0.204257, 0.22407, 0.243493, 0.264567, 0.284231, 
0.305546, 0.326471, 0.345985, 0.367151, 0.386907, 0.408314, 0.429331, 
0.448938, 0.470196, 0.490044, 0.509502, 0.530611, 0.55031, 0.57166, 
0.59262, 0.61217, 0.633371, 0.653162, 0.672563, 0.693616, 0.713258, 
0.734551, 0.754434, 0.773927, 0.795071, 0.814805, 0.83619, 0.857186, 
0.876771, 0.898007, 0.917833, 0.93727, 0.958357, 0.978034, 0.978377, 
0.978721, 0.979407, 0.98078, 0.983526, 0.989017, 0.98936, 0.989704, 
0.99039, 0.991763, 0.994509, 0.994852, 0.995195, 0.995881, 0.997254, 
0.997597, 0.997941, 0.998627, 0.99897, 0.999314, 0.999657, 1.],
y:[1.5708, 1.5708, 1.5708, 1.57079, 1.57079, 1.57079, 1.57079, 
1.57079, 1.57078, 1.57078, 1.57078, 1.57077, 1.57077, 1.57077, 
1.57076, 1.57075, 1.57074, 1.57074, 1.57072, 1.57072, 1.57071, 
1.57069, 1.57068, 1.57067, 1.57065, 1.57061, 1.5706, 1.57059, 
1.57056, 1.5705, 1.57049, 1.57047, 1.57044, 1.57038, 1.57036, 
1.57034, 1.5703, 1.57022, 1.57005, 1.57002, 1.57, 1.56995, 1.56983, 
1.56959, 1.56956, 1.56953, 1.56946, 1.56932, 1.56902, 1.56898, 
1.56894, 1.56886, 1.5687, 1.56834, 1.56755, 1.56751, 1.56746, 
1.56735, 1.56715, 1.56672, 1.56579, 1.56573, 1.56567, 1.56554, 
1.56529, 1.56476, 1.56364, 1.56357, 1.5635, 1.56335, 1.56306, 
1.56245, 1.56117, 1.56109, 1.561, 1.56083, 1.56049, 1.5598, 1.55833, 
1.55823, 1.55812, 1.55792, 1.5575, 1.55664, 1.55484, 1.55092, 1.5508, 
1.55068, 1.55043, 1.54994, 1.54894, 1.54687, 1.54246, 1.53162, 
1.5193, 1.50632, 1.49064, 1.4745, 1.45729, 1.43709, 1.41682, 1.39332, 
1.36874, 1.34447, 1.31673, 1.28952, 1.25863, 1.22693, 1.19614, 
1.16148, 1.12795, 1.09402, 1.05605, 1.01957, 0.978929, 0.93796, 
0.89883, 0.855438, 0.814076, 0.772764, 0.727126, 0.683827, 0.636155, 
0.590996, 0.546162, 0.496954, 0.450532, 0.399738, 0.349432, 0.30216, 
0.250582, 0.202176, 0.154531, 0.102677, 0.0541875, 0.0533412, 
0.0524948, 0.050802, 0.0474163, 0.0406442, 0.0270978, 0.0262511, 
0.0254044, 0.0237109, 0.0203238, 0.0135494, 0.0127026, 0.0118558, 
0.0101621, 0.00677481, 0.00592797, 0.00508112, 0.00338744, 
0.00254059, 0.00169374, 0.000846897, 5.03551e-8],
mode: 'lines',
name: 'standard'

};

var plot1 = {
x:[2.04082e-8, 0.000306718, 0.000613415, 0.000920113, 0.00122681, 
0.00153351, 0.00184021, 0.0021469, 0.0024536, 0.0027603, 0.003067, 
0.00337369, 0.00368039, 0.00398709, 0.00429379, 0.00490718, 
0.00521388, 0.00552058, 0.00613397, 0.00644067, 0.00674737, 
0.00736076, 0.00766746, 0.00797416, 0.00858755, 0.00981434, 0.010121, 
0.0104277, 0.0110411, 0.0122679, 0.0125746, 0.0128813, 0.0134947, 
0.0147215, 0.0150282, 0.0153349, 0.0159483, 0.0171751, 0.0196287, 
0.0199612, 0.0202937, 0.0209586, 0.0222886, 0.0249486, 0.0252811, 
0.0256136, 0.0262786, 0.0276085, 0.0302685, 0.030601, 0.0309335, 
0.0315985, 0.0329285, 0.0355884, 0.0409084, 0.0412188, 0.0415293, 
0.0421502, 0.043392, 0.0458757, 0.0508431, 0.0511536, 0.051464, 
0.052085, 0.0533268, 0.0558105, 0.0607779, 0.0610823, 0.0613866, 
0.0619954, 0.0632129, 0.0656478, 0.0705178, 0.0708221, 0.0711265, 
0.0717352, 0.0729527, 0.0753877, 0.0802576, 0.0805878, 0.080918, 
0.0815783, 0.082899, 0.0855404, 0.0908231, 0.101388, 0.101697, 
0.102005, 0.102621, 0.103854, 0.106319, 0.111249, 0.121109, 0.142481, 
0.163463, 0.183035, 0.204257, 0.22407, 0.243493, 0.264567, 0.284231, 
0.305546, 0.326471, 0.345985, 0.367151, 0.386907, 0.408314, 0.429331, 
0.448938, 0.470196, 0.490044, 0.509502, 0.530611, 0.55031, 0.57166, 
0.59262, 0.61217, 0.633371, 0.653162, 0.672563, 0.693616, 0.713258, 
0.734551, 0.754434, 0.773927, 0.795071, 0.814805, 0.83619, 0.857186, 
0.876771, 0.898007, 0.917833, 0.93727, 0.958357, 0.978034, 0.978377, 
0.978721, 0.979407, 0.98078, 0.983526, 0.989017, 0.98936, 0.989704, 
0.99039, 0.991763, 0.994509, 0.994852, 0.995195, 0.995881, 0.997254, 
0.997597, 0.997941, 0.998627, 0.99897, 0.999314, 0.999657, 1.],
y:[1.11072, 1.11072, 1.11072, 1.11072, 1.11072, 1.11072, 1.11072, 
1.11072, 1.11072, 1.11072, 1.11072, 1.11072, 1.11072, 1.11072, 
1.11071, 1.11071, 1.11071, 1.11071, 1.11071, 1.11071, 1.11071, 
1.1107, 1.1107, 1.1107, 1.1107, 1.11069, 1.11069, 1.11068, 1.11068, 
1.11067, 1.11067, 1.11066, 1.11066, 1.11065, 1.11064, 1.11064, 
1.11063, 1.11062, 1.11059, 1.11058, 1.11058, 1.11057, 1.11055, 
1.11051, 1.1105, 1.1105, 1.11048, 1.11046, 1.11041, 1.1104, 1.11039, 
1.11038, 1.11035, 1.11029, 1.11015, 1.11014, 1.11013, 1.11011, 
1.11008, 1.11, 1.10984, 1.10982, 1.10981, 1.10979, 1.10975, 1.10965, 
1.10946, 1.10944, 1.10943, 1.1094, 1.10935, 1.10924, 1.10902, 1.109, 
1.10899, 1.10896, 1.1089, 1.10877, 1.10851, 1.1085, 1.10848, 1.10844, 
1.10837, 1.10822, 1.1079, 1.1072, 1.10718, 1.10716, 1.10712, 1.10703, 
1.10685, 1.10648, 1.1057, 1.10377, 1.10158, 1.09926, 1.09646, 
1.09357, 1.09047, 1.08683, 1.08316, 1.07889, 1.07441, 1.06996, 
1.06486, 1.05983, 1.05409, 1.04817, 1.04239, 1.03584, 1.02946, 
1.02297, 1.01566, 1.00858, 1.00064, 0.992566, 0.984794, 0.976105, 
0.967748, 0.95933, 0.949943, 0.940951, 0.93095, 0.921376, 0.911772, 
0.901112, 0.89094, 0.879674, 0.868373, 0.857618, 0.845727, 0.834413, 
0.823125, 0.810661, 0.79883, 0.798622, 0.798414, 0.797998, 0.797164, 
0.795494, 0.792144, 0.791934, 0.791724, 0.791304, 0.790463, 0.788778, 
0.788567, 0.788357, 0.787935, 0.78709, 0.786879, 0.786667, 0.786245, 
0.786033, 0.785821, 0.78561, 0.785398],
mode: 'lines',
name: 'half angle'
};

var plot2 = {
x:[2.04082e-8, 0.000306718, 0.000613415, 0.000920113, 0.00122681, 
0.00153351, 0.00184021, 0.0021469, 0.0024536, 0.0027603, 0.003067, 
0.00337369, 0.00368039, 0.00398709, 0.00429379, 0.00490718, 
0.00521388, 0.00552058, 0.00613397, 0.00644067, 0.00674737, 
0.00736076, 0.00766746, 0.00797416, 0.00858755, 0.00981434, 0.010121, 
0.0104277, 0.0110411, 0.0122679, 0.0125746, 0.0128813, 0.0134947, 
0.0147215, 0.0150282, 0.0153349, 0.0159483, 0.0171751, 0.0196287, 
0.0199612, 0.0202937, 0.0209586, 0.0222886, 0.0249486, 0.0252811, 
0.0256136, 0.0262786, 0.0276085, 0.0302685, 0.030601, 0.0309335, 
0.0315985, 0.0329285, 0.0355884, 0.0409084, 0.0412188, 0.0415293, 
0.0421502, 0.043392, 0.0458757, 0.0508431, 0.0511536, 0.051464, 
0.052085, 0.0533268, 0.0558105, 0.0607779, 0.0610823, 0.0613866, 
0.0619954, 0.0632129, 0.0656478, 0.0705178, 0.0708221, 0.0711265, 
0.0717352, 0.0729527, 0.0753877, 0.0802576, 0.0805878, 0.080918, 
0.0815783, 0.082899, 0.0855404, 0.0908231, 0.101388, 0.101697, 
0.102005, 0.102621, 0.103854, 0.106319, 0.111249, 0.121109, 0.142481, 
0.163463, 0.183035, 0.204257, 0.22407, 0.243493, 0.264567, 0.284231, 
0.305546, 0.326471, 0.345985, 0.367151, 0.386907, 0.408314, 0.429331, 
0.448938, 0.470196, 0.490044, 0.509502, 0.530611, 0.55031, 0.57166, 
0.59262, 0.61217, 0.633371, 0.653162, 0.672563, 0.693616, 0.713258, 
0.734551, 0.754434, 0.773927, 0.795071, 0.814805, 0.83619, 0.857186, 
0.876771, 0.898007, 0.917833, 0.93727, 0.958357, 0.978034, 0.978377, 
0.978721, 0.979407, 0.98078, 0.983526, 0.989017, 0.98936, 0.989704, 
0.99039, 0.991763, 0.994509, 0.994852, 0.995195, 0.995881, 0.997254, 
0.997597, 0.997941, 0.998627, 0.99897, 0.999314, 0.999657, 1.],
y:[0.785398, 0.785398, 0.785398, 0.785399, 0.785399, 0.785399, 
0.7854, 0.7854, 0.785401, 0.785402, 0.785403, 0.785404, 0.785405, 
0.785406, 0.785407, 0.78541, 0.785411, 0.785413, 0.785416, 0.785418, 
0.78542, 0.785424, 0.785427, 0.785429, 0.785434, 0.785445, 0.785448, 
0.785451, 0.785457, 0.785471, 0.785475, 0.785479, 0.785486, 0.785503, 
0.785508, 0.785512, 0.785521, 0.785541, 0.785585, 0.785591, 0.785598, 
0.785611, 0.785639, 0.7857, 0.785708, 0.785716, 0.785733, 0.785768, 
0.785842, 0.785852, 0.785862, 0.785882, 0.785924, 0.786012, 0.786209, 
0.786222, 0.786234, 0.78626, 0.786311, 0.786419, 0.786652, 0.786667, 
0.786683, 0.786714, 0.786777, 0.786909, 0.787191, 0.787209, 0.787227, 
0.787263, 0.787337, 0.78749, 0.787812, 0.787833, 0.787854, 0.787897, 
0.787982, 0.788158, 0.788527, 0.788553, 0.788579, 0.788631, 0.788737, 
0.788954, 0.789408, 0.790399, 0.79043, 0.790461, 0.790522, 0.790647, 
0.7909, 0.791425, 0.792547, 0.795316, 0.798487, 0.801855, 0.805963, 
0.810233, 0.814837, 0.82031, 0.825875, 0.832422, 0.839385, 0.846372, 
0.854504, 0.862632, 0.872048, 0.881931, 0.891744, 0.903056, 0.914275, 
0.925918, 0.939302, 0.95253, 0.967712, 0.983516, 0.999099, 1.01697, 
1.03461, 1.05284, 1.07375, 1.09437, 1.11801, 1.14137, 1.16554, 
1.19328, 1.22068, 1.25215, 1.28495, 1.31738, 1.35469, 1.39168, 
1.43011, 1.47442, 1.51842, 1.51921, 1.52, 1.52159, 1.52477, 1.53118, 
1.54416, 1.54498, 1.5458, 1.54744, 1.55073, 1.55736, 1.5582, 1.55903, 
1.5607, 1.56405, 1.56489, 1.56573, 1.56742, 1.56826, 1.5691, 1.56995, 
1.5708],
mode: 'lines',
name: 'Cayley'
};

var plot3 = { x:[2.04082e-8, 1.], y:[1., 1.], mode: 'lines', name: 'log' };

var plot4 = {
x:[2.04082e-8, 0.000306718, 0.000613415, 0.000920113, 0.00122681, 
0.00153351, 0.00184021, 0.0021469, 0.0024536, 0.0027603, 0.003067, 
0.00337369, 0.00368039, 0.00398709, 0.00429379, 0.00490718, 
0.00521388, 0.00552058, 0.00613397, 0.00644067, 0.00674737, 
0.00736076, 0.00766746, 0.00797416, 0.00858755, 0.00981434, 0.010121, 
0.0104277, 0.0110411, 0.0122679, 0.0125746, 0.0128813, 0.0134947, 
0.0147215, 0.0150282, 0.0153349, 0.0159483, 0.0171751, 0.0196287, 
0.0199612, 0.0202937, 0.0209586, 0.0222886, 0.0249486, 0.0252811, 
0.0256136, 0.0262786, 0.0276085, 0.0302685, 0.030601, 0.0309335, 
0.0315985, 0.0329285, 0.0355884, 0.0409084, 0.0412188, 0.0415293, 
0.0421502, 0.043392, 0.0458757, 0.0508431, 0.0511536, 0.051464, 
0.052085, 0.0533268, 0.0558105, 0.0607779, 0.0610823, 0.0613866, 
0.0619954, 0.0632129, 0.0656478, 0.0705178, 0.0708221, 0.0711265, 
0.0717352, 0.0729527, 0.0753877, 0.0802576, 0.0805878, 0.080918, 
0.0815783, 0.082899, 0.0855404, 0.0908231, 0.101388, 0.101697, 
0.102005, 0.102621, 0.103854, 0.106319, 0.111249, 0.121109, 0.142481, 
0.163463, 0.183035, 0.204257, 0.22407, 0.243493, 0.264567, 0.284231, 
0.305546, 0.326471, 0.345985, 0.367151, 0.386907, 0.408314, 0.429331, 
0.448938, 0.470196, 0.490044, 0.509502, 0.530611, 0.55031, 0.57166, 
0.59262, 0.61217, 0.633371, 0.653162, 0.672563, 0.693616, 0.713258, 
0.734551, 0.754434, 0.773927, 0.795071, 0.814805, 0.83619, 0.857186, 
0.876771, 0.898007, 0.917833, 0.93727, 0.958357, 0.978034, 0.978377, 
0.978721, 0.979407, 0.98078, 0.983526, 0.989017, 0.98936, 0.989704, 
0.99039, 0.991763, 0.994509, 0.994852, 0.995195, 0.995881, 0.997254, 
0.997597, 0.997941, 0.998627, 0.99897, 0.999314, 0.999657, 1.],
y:[1.02617, 1.02617, 1.02617, 1.02617, 1.02617, 1.02617, 1.02617, 
1.02617, 1.02617, 1.02617, 1.02617, 1.02617, 1.02617, 1.02617, 
1.02617, 1.02617, 1.02617, 1.02617, 1.02617, 1.02617, 1.02617, 
1.02617, 1.02617, 1.02617, 1.02617, 1.02616, 1.02616, 1.02616, 
1.02616, 1.02616, 1.02616, 1.02616, 1.02616, 1.02616, 1.02615, 
1.02615, 1.02615, 1.02615, 1.02614, 1.02614, 1.02614, 1.02614, 
1.02613, 1.02612, 1.02612, 1.02612, 1.02612, 1.02611, 1.0261, 1.0261, 
1.0261, 1.02609, 1.02609, 1.02607, 1.02604, 1.02604, 1.02604, 
1.02603, 1.02602, 1.02601, 1.02597, 1.02597, 1.02596, 1.02596, 
1.02595, 1.02593, 1.02588, 1.02588, 1.02587, 1.02587, 1.02586, 
1.02583, 1.02578, 1.02578, 1.02577, 1.02577, 1.02575, 1.02572, 
1.02566, 1.02566, 1.02565, 1.02565, 1.02563, 1.02559, 1.02552, 
1.02536, 1.02535, 1.02535, 1.02534, 1.02532, 1.02528, 1.02519, 
1.02501, 1.02457, 1.02406, 1.02352, 1.02287, 1.0222, 1.02148, 
1.02064, 1.01979, 1.01879, 1.01775, 1.01672, 1.01552, 1.01435, 
1.01301, 1.01162, 1.01027, 1.00873, 1.00723, 1.0057, 1.00398, 1.0023, 
1.00042, 0.998509, 0.996663, 0.994594, 0.992601, 0.990589, 0.98834, 
0.986181, 0.983775, 0.981465, 0.979143, 0.97656, 0.974088, 0.971343, 
0.968581, 0.965946, 0.963023, 0.960235, 0.957445, 0.954354, 0.951412, 
0.95136, 0.951308, 0.951204, 0.950996, 0.95058, 0.949744, 0.949692, 
0.94964, 0.949535, 0.949325, 0.948904, 0.948851, 0.948799, 0.948693, 
0.948482, 0.94843, 0.948377, 0.948271, 0.948218, 0.948165, 0.948112,
0.948059],
mode: 'lines',
name: 'half angle/half angle'
};

var plot5 = {
x:[2.04082e-8, 0.000306718, 0.000613415, 0.000920113, 0.00122681, 
0.00153351, 0.00184021, 0.0021469, 0.0024536, 0.0027603, 0.003067, 
0.00337369, 0.00368039, 0.00398709, 0.00429379, 0.00490718, 
0.00521388, 0.00552058, 0.00613397, 0.00644067, 0.00674737, 
0.00736076, 0.00766746, 0.00797416, 0.00858755, 0.00981434, 0.010121, 
0.0104277, 0.0110411, 0.0122679, 0.0125746, 0.0128813, 0.0134947, 
0.0147215, 0.0150282, 0.0153349, 0.0159483, 0.0171751, 0.0196287, 
0.0199612, 0.0202937, 0.0209586, 0.0222886, 0.0249486, 0.0252811, 
0.0256136, 0.0262786, 0.0276085, 0.0302685, 0.030601, 0.0309335, 
0.0315985, 0.0329285, 0.0355884, 0.0409084, 0.0412188, 0.0415293, 
0.0421502, 0.043392, 0.0458757, 0.0508431, 0.0511536, 0.051464, 
0.052085, 0.0533268, 0.0558105, 0.0607779, 0.0610823, 0.0613866, 
0.0619954, 0.0632129, 0.0656478, 0.0705178, 0.0708221, 0.0711265, 
0.0717352, 0.0729527, 0.0753877, 0.0802576, 0.0805878, 0.080918, 
0.0815783, 0.082899, 0.0855404, 0.0908231, 0.101388, 0.101697, 
0.102005, 0.102621, 0.103854, 0.106319, 0.111249, 0.121109, 0.142481, 
0.163463, 0.183035, 0.204257, 0.22407, 0.243493, 0.264567, 0.284231, 
0.305546, 0.326471, 0.345985, 0.367151, 0.386907, 0.408314, 0.429331, 
0.448938, 0.470196, 0.490044, 0.509502, 0.530611, 0.55031, 0.57166, 
0.59262, 0.61217, 0.633371, 0.653162, 0.672563, 0.693616, 0.713258, 
0.734551, 0.754434, 0.773927, 0.795071, 0.814805, 0.83619, 0.857186, 
0.876771, 0.898007, 0.917833, 0.93727, 0.958357, 0.978034, 0.978377, 
0.978721, 0.979407, 0.98078, 0.983526, 0.989017, 0.98936, 0.989704, 
0.99039, 0.991763, 0.994509, 0.994852, 0.995195, 0.995881, 0.997254, 
0.997597, 0.997941, 0.998627, 0.99897, 0.999314, 0.999657, 1.],
y:[0.948059, 0.948059, 0.94806, 0.94806, 0.94806, 0.94806, 0.94806, 
0.94806, 0.94806, 0.948061, 0.948061, 0.948061, 0.948061, 0.948062, 
0.948062, 0.948063, 0.948063, 0.948064, 0.948065, 0.948066, 0.948066, 
0.948067, 0.948068, 0.948069, 0.94807, 0.948074, 0.948074, 0.948075, 
0.948077, 0.948081, 0.948083, 0.948084, 0.948086, 0.948091, 0.948092, 
0.948094, 0.948097, 0.948103, 0.948116, 0.948118, 0.94812, 0.948124, 
0.948132, 0.94815, 0.948153, 0.948155, 0.94816, 0.948171, 0.948193, 
0.948196, 0.948199, 0.948205, 0.948218, 0.948245, 0.948304, 0.948308, 
0.948312, 0.948319, 0.948335, 0.948367, 0.948437, 0.948442, 0.948447, 
0.948456, 0.948475, 0.948515, 0.9486, 0.948605, 0.948611, 0.948622, 
0.948644, 0.94869, 0.948787, 0.948793, 0.948799, 0.948812, 0.948838, 
0.948891, 0.949002, 0.94901, 0.949017, 0.949033, 0.949065, 0.94913, 
0.949266, 0.949564, 0.949573, 0.949582, 0.949601, 0.949638, 0.949714, 
0.949871, 0.950207, 0.951034, 0.951977, 0.952974, 0.954185, 0.955438, 
0.956781, 0.958367, 0.95997, 0.961841, 0.963815, 0.965778, 0.968044, 
0.970287, 0.972858, 0.975527, 0.978147, 0.981132, 0.984054, 0.987049, 
0.990444, 0.993752, 0.99749, 1.00132, 1.00503, 1.00922, 1.01327, 
1.01739, 1.02203, 1.02651, 1.03154, 1.0364, 1.04132, 1.04684, 
1.05216, 1.05811, 1.06415, 1.06997, 1.07647, 1.08273, 1.08905, 
1.09612, 1.10291, 1.10303, 1.10315, 1.10339, 1.10387, 1.10484, 
1.10678, 1.10691, 1.10703, 1.10727, 1.10776, 1.10874, 1.10887, 
1.10899, 1.10924, 1.10973, 1.10985, 1.10998, 1.11023, 1.11035, 
1.11047, 1.1106, 1.11072],
mode: 'lines',
name: 'half angle/Cayley'
};

var data = [plot0, plot1,plot2,plot3,plot4,plot5];

var layout = {
  title:  'density function',
  xaxis: { nticks: 10 },
  yaxis: { nticks: 20 },
  height: 376,
  width:  626
};

Plotly.newPlot('fig3', data, layout, {displaylogo: false, autosizable: true});
</script>


