---
layout:       post
title:        Volume element of SO(3) and average uniform random rotation angle
categories:   [quaternions]
tags:         [rotation, quantization, distribution, sampling, map]
description:  Find the expression of a volume element of rotations.
plotly: true
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

------

Attempt at human readable summary
------

\\
I'm going to try summarize this post for a wide target audience.  Let's say we wanted to find the average angle of a uniform random rotation in two dimensions.  This is the same as points uniformly distributed on the surface (arc lengths) of the unit circle and if we measure angles as minimum magnitude then the answer is the probably unsurprising value of  $\frac{\pi}{2}$ radians.  In three dimensions this is the same problem as uniform random rotations about a fixed axis.  If we extend this to allow any axis of rotation bound to some plane (say the $\lbrace xy\rbrace$-plane) then the result is the same.  This is equivalent to points uniformly distributed on the surface of the unit 3D sphere, choose a reference direction (logical choice is $+\mathbf{z}$) and finding the rotation that maps that direction $\left(\mathbf{z}\right)$ to the random point.  I've skimmed the math of random generation a previous post [Uniform points on disc, circle, sphere and caps]({{site.base}}/distribution/2016/11/28/Uniform.html).  Additionally the median of both of these is $\frac{\pi}{2}$ as well.  To be clear on terminology the average is: "let's generate $n$ values and average their result and let $n$ go to infinity" and median is "let's generate $n$ results, sort them and find the value at $\frac{n}{2}$".  So in terms of uniform distributed point-sets the median is the value (here angle) where the length/area/volume is divided in two.


Moving to the general case of 3D uniform random rotations is equivalent to points uniformly distributed on the 4D unit hemisphere.  This is hard to visualize so let's look at a 3D slice of it:

![stretch]({{site.base}}/assets/figures/misc/VolSo3.png '4D half sphere'){: .center-image }

\\
The vertical axis is line of reals (aka scalars, typically called 'w') and the "north pole" is the quaternion value of one which is the identity (no rotation) transform.  The other two dimensions can be any planar slice of what we'd typically think of as normal 3D space, but let's say that it's the $\lbrace xy\rbrace$-plane. This choice means that hemisphere we're visualizing is the set of all rotations with an axis in the $\lbrace xy\rbrace$-plane, so the second example from above.  This specialized case actually is a 3D hemisphere (the $z$ component of the quaternion is always zero) and the "white-ish" circle indicates the $\frac{\pi}{2}$ split which is at an angle of $\frac{\pi}{4}$ with respect to the real (vertical) axis because of the half/double angle relationship between 3D rotation angle and quaternion space measurements.

<div id="cap" style="width:100%"></div>


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
Great it agrees with what we could have looked up on Wikipedia (EDIT: [MathWorld](http://mathworld.wolfram.com/Hypersphere.html) has the volume element)


<br>

------

Volume of $SO\left(3\right)$ <small>take 2</small>
------

\\
I don't need it but since I'm typing we could have done this a second way. ([Wikipedia](https://en.wikipedia.org/wiki/Volume_element#Volume_element_of_manifolds), [Encyclopedia of Mathematics](https://www.encyclopediaofmath.org/index.php/Jacobian))  Without adding an extra dimension the Jacobian is:

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
then the [metric tensor](http://en.wikipedia.org/wiki/Metric_tensor):

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
Not explictly mentioned so far is that $\theta$ is the measure in quaternion space which translates into a $2\theta$ rotation. I want to find the average rotation angle $\phi$ of a uniform distribution of random rotations.  We have $q_w = \cos\left(\frac{\phi}{2}\right) = \cos\left(\theta\right)$ and the [expected value](http://en.wikipedia.org/wiki/Expected_value) is:

$$ \begin{eqnarray*}
E\lbrack\phi\rbrack & = & 2~E\lbrack\theta\rbrack     \\
                    & = & \frac{1}{A}\int 2~\theta~dA \\
                    & = & \frac{2}{\pi^2}\int ~\theta~\sin\left(\theta\right)^2~\sin\left(\alpha\right)~d\theta~d\alpha~d\beta \\
                    & = & \frac{4}{\pi^2} \left(\pi + \frac{\pi^3}{4}\right) \\
                    & = & \frac{\pi}{2} + \frac{2}{\pi}
\end{eqnarray*} $$

There's surely an easier way to do this and I'm shaky on the math here, so an empirical validation: [CLICK](https://gist.github.com/Marc-B-Reynolds/68ad708c950f57f0e38a445f9e9ef697)

EDIT: Same result discussed in this [math.stackexchange](https://math.stackexchange.com/questions/464419/mean-value-of-the-rotation-angle-is-126-5%C2%B0) post based on this paper ([PDF](http://www.theworld.com/~sweetser/quaternions/ps/stanfordaiwp79-salamin.pdf)).  Also I was using probablity verbage here where I should have used geometric.  I'm finding the "center-of-mass" of the point set since the "average" is summing up the angles and dividing by number of points (extended from discrete to continous). 

<br>

------

Cumulative volume of angle range <small>volume of cap & median random rotation<small>
------

\\
Finding the volume of all 3D rotations with angle $\phi$ and less is simply changing the upper bound of the $\theta$ integral:

$$ \int_{0}^{2\pi} \int_{0}^{\pi} \int_{0}^{\phi/2} \sin\left(\theta\right)^2~\sin\left(\alpha\right)~d\theta~d\alpha~d\beta = \pi\left(\phi - \sin\left(\phi\right)\right) $$

\\
Solving the above for half the total volume of rotations gives $ \phi \approx 2.30988 $ radians (~132.35 degrees) which is the median angle of uniform random rotations.

<br>


<script>

'use strict';

const NUM_SAMPLES = 100;
const IPI = 1.0/Math.PI;

var x_coord = new Array(NUM_SAMPLES);
var y_cap   = new Array(NUM_SAMPLES);

{
  let K = Math.PI/(NUM_SAMPLES-1.0);

  for(var i=0; i<NUM_SAMPLES; i++) {
    let t = K*i
    x_coord[i] = t;
	y_cap[i]   = Math.PI*(t-Math.sin(t));
  }
}

var median    = {y:[4.934802200544], x:[2.309881460010058], name:'median', mode:'markers'};
var average   = {y:[4.408616671502], x:[2.207416099162478], name:'average', mode:'markers'};
var cap_trace = {y:y_cap, x:x_coord, name:'cumulative volume'};
var cap_data  = [cap_trace,median,average];



var cap_layout = {
  title:  'cumulative volume',
  xaxis: { nticks: 10 },
  yaxis: { nticks: 12 },
  height: 400,
  width:  800
};

Plotly.newPlot('cap', cap_data, cap_layout, {displaylogo: false, autosizable: true});

</script>
