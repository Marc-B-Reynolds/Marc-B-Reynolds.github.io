---
layout: post
title:  Converting to Euler & Tait-Bryan
categories: [math]
tags : [euler, rotations]
description : A note on reducing errors when converting to Euler angles or Tait-Bryan paramerters.
plotly:       true
---

\\
In the process of working on a quaternion quantization post I needed to round-trip between the semi-standard yaw/pitch/roll representation and quaternions. In the process I discovered just how inaccurate my existing conversion code and all those that I could find on the web are.  This is note about the duct-tape & super-glue I applied to lower the round-trip errors to an acceptable range for my purposes.  I'm not going to carry through beyond that point.  I set myself a target peak error of ~0.00015 degrees.  I'll focus on converting from quaternions but will make some nods to matrices.

Aside: there are a number of methods to map a quaternion to three values which are computationally cheaper and introduce only small errors to round-trip.

Test framework can be found: [HERE](http://github.com/Marc-B-Reynolds/Stand-alone-junk/blob/master/src/Posts/tait2q.c)
<br>

------

Prelim
------

\\
Informally *Euler angles* refer to a parameterization of rotations which is a set of three angles about three predefined axes. It's usually noted that this gives us 12 possible choices $\left(3\times2\times2\right)$, three for the first axis and two choices for the second and third (since we can't repeat the immediately previous choice).  And then that number is doubled for left vs. right handed systems.  And we could double that again for fixed frame (global rotation view/space-fixed) or inertial frame (local rotation view/body-fixed).  I'd argue that there are really only two distinct versions: All three axes choices are present (Tait-Bryan) and the first and last are repeats (proper Euler angles).  Choice of axis/handedness is just sign flips/permutations and local vs. global is the choice of reading the componsition left-to-right vs. right-to-left.


<br>

------

Tait-Bryan: ZYX <small></small>
------

Seemly the most common choice is that of Tait-Bryan ZYX.  This wikipedia picture illustrates the parameterization (except negate the y-axis and reverse the direction of $\phi$...use your imagination):


![rot](https://upload.wikimedia.org/wikipedia/commons/3/30/Plane_with_ENU_embedded_axes.svg 'there you go'){: width="50%" .center-image }


First I'll walk through the basics from a quaternion perspective.  Define a standard basis set:

$$
\begin{eqnarray*}
  \mathbf{x} & = & \left(1,0,0\right) \\ 
  \mathbf{y} & = & \left(0,1,0\right) \\ 
  \mathbf{z} & = & \left(0,0,1\right) 
\end{eqnarray*}
$$

We take $\mathbf{z}$ to be up, the default "facing" is $\mathbf{x}$ and (of course) right-handed, then the specific parameterization is:

1.  rotate about local $\mathbf{z}$ by $\psi$ (*yaw*)

2.  rotate about local $\mathbf{y}$ by $\theta$ (*pitch*)

3.  rotate about local $\mathbf{x}$ by $\phi$ (*roll*)

\\
Then the quaternion expression of interest is:

$$ \begin{equation} \label{eq:tait}
Q = w + \left(x,~y,~z\right) = e^{\frac{\psi}{2}\mathbf{z}}e^{\frac{\theta}{2}\mathbf{y}}e^{\frac{\phi}{2}\mathbf{x}}
\end{equation} $$

and recall that:

$$
  e^{\alpha\mathbf{u}} = \cos\left(\alpha\right) + \sin\left(\alpha\right)\mathbf{u}
$$

<br>

------

ZYX to quaternion <small>and matrix</small>
------

\\
To formulate the Euler to quaternion conversion we simply expand the right most expression of $\eqref{eq:tait}$ and set to the middle expression:

$$
\begin{eqnarray}
  w &=& \cos\left(\frac{\phi}{2}\right) \cos\left(\frac{\theta}{2}\right) \cos\left(\frac{\psi}{2}\right) + \sin\left(\frac{\phi}{2}\right) \sin\left(\frac{\theta}{2}\right) \sin\left(\frac{\psi}{2}\right) \label{qw} \\
  x &=& \sin\left(\frac{\phi}{2}\right) \cos\left(\frac{\theta}{2}\right) \cos\left(\frac{\psi}{2}\right) - \cos\left(\frac{\phi}{2}\right) \sin\left(\frac{\theta}{2}\right) \sin\left(\frac{\psi}{2}\right) \label{qx} \\ 
  y &=& \cos\left(\frac{\phi}{2}\right) \sin\left(\frac{\theta}{2}\right) \cos\left(\frac{\psi}{2}\right) + \sin\left(\frac{\phi}{2}\right) \cos\left(\frac{\theta}{2}\right) \sin\left(\frac{\psi}{2}\right) \label{qy} \\ 
  z &=&  \cos\left(\frac{\phi}{2}\right) \cos\left(\frac{\theta}{2}\right) \sin\left(\frac{\psi}{2}\right)-\sin\left(\frac{\phi}{2}\right) \sin\left(\frac{\theta}{2}\right) \cos\left(\frac{\psi}{2}\right) \label{qz}

\end{eqnarray}
$$

\\
(really the matrix conversion is hidden in the next section.)

<br>

------

And back again... <small>as if we can compute with Reals</small>
------

\\
Given $Q$ we can use the similarity transform to find $\left(\mathbf{x}',~\mathbf{y}',~\mathbf{z}'\right)$: 

$$
\begin{eqnarray}
  \mathbf{x}' &=& Q\mathbf{x}Q^* \nonumber \\
              &=& \left(w^2+x^2-y^2-z^2, ~2xy+2wz, ~2xz-2wy\right) \label{xpu} \\
              &=& \left(1-2\left(y^2+z^2\right), ~2(xy+wz), ~2(xz-wy)\right) \label{xp} \\
              &=& \left(\cos(\theta) \cos(\psi), ~\cos(\theta) \sin(\psi),~ -\sin(\theta)\right) \label{xpt} \\
              &=& \left(m_{00}, ~m_{01}, ~m_{02} \right) \nonumber \\
			  \nonumber \\
  \mathbf{y}' &=& Q\mathbf{y}Q^* \nonumber \\
              &=& \left(2xy-2wz, ~w^2-x^2+y^2-z^2, ~2wx+2yz\right) \label{ypu} \\
              &=& \left(2(xy-wz), ~1-2\left(x^2+z^2\right), ~2(wx+yz)\right) \label{yp} \\
 &=& \left(\sin(\theta)\cos(\psi)\sin(\phi)-\sin(\psi)\cos(\phi),~\sin(\theta)\sin(\psi)\sin(\phi)+\cos(\psi)\cos(\phi),~\cos(\theta)\sin(\phi)\right) \label{ypt} \\
              &=& \left(m_{10}, ~m_{11}, ~m_{12} \right) \nonumber \\
			  \nonumber \\
  \mathbf{z}' &=& Q\mathbf{z}Q^* \nonumber \\
              &=& \left(2wy+2xz, ~2yz-2wx, ~w^2-x^2-y^2+z^2\right) \label{zpu} \\
              &=& \left(2(wy+xz), ~2(yz-wx), ~1-2\left(x^2+y^2\right)\right) \label{zp} \\
 &=& \left(\sin(\theta)\cos(\psi)\cos(\phi)+\sin(\psi)\sin(\phi),~\sin(\theta)\sin(\psi)\cos(\phi)-\cos(\psi)\sin(\phi),~\cos(\theta)\cos(\phi)\right) \label{zpt} \\
              &=& \left(m_{20}, ~m_{21}, ~m_{22} \right) \nonumber \\
\end{eqnarray}
$$

\\
Where the $m_{ij}$ are to show the equivalent matrix $R$.

$$
R = 
\left(
\begin{matrix}
  \mathbf{x}' \\
  \mathbf{y}' \\
  \mathbf{z}' \\
\end{matrix}
\right) =
\left(
\begin{matrix}
  m_{00} & m_{01} & m_{02} \\
  m_{10} & m_{11} & m_{12} \\
  m_{20} & m_{21} & m_{22} \\
\end{matrix}
\right)
$$


Ignoring the degenerate region we can easily find $\theta$ and $\psi$ from $ \mathbf{x}'$:

$$
\begin{eqnarray}
  \theta &=& \text{asin} \left(-m_{02}\right) \label{thetam} \\
         &=& \text{asin} \left(-2\left(xz-wy\right)\right) \label{thetaq} \\
  \nonumber \\
  \psi &=& \text{atan2}\left(m_{01},~m_{00} \right)  \label{psim} \\
	   &=& \text{atan2}\left(xy+wz,~\frac{1}{2}-\left(y^2+z^2\right)\right)      \label{psiq2}
\end{eqnarray}
$$

\\
and $\phi$ from the z components of $\mathbf{y}'$ and $\mathbf{z}'$

$$
\begin{eqnarray}
  \phi  &=&  \text{atan2}\left(m_{12},~m_{22}\right) \\
        &=&  \text{atan2}\left(yz+wx,~\frac{1}{2}-\left(x^2+y^2\right)\right) \label{phiq1}
\end{eqnarray}
$$

\\
The above computations explode into an infinite number of solutions when $\cos\left(\theta\right)=0$
($\theta=\pm\frac{\pi}{2}$).  The matrix form then becomes:

$$
\begin{eqnarray*}
\theta &=& \frac{\pi}{2} &\to& \left(
\begin{matrix}
  0 & 0 & -1 \\
  \sin\left( \phi - \psi \right) &  \cos\left( \phi - \psi \right) & 0 \\
  \cos\left( \phi - \psi \right) & -\sin\left( \phi - \psi \right) & 0 \\
\end{matrix}
\right) \\
\theta &=& -\frac{\pi}{2} &\to& \left(
\begin{matrix}
  0 & 0 & 1 \\
  -\sin\left(\phi + \psi \right) &  \cos\left( \phi + \psi \right) & 0 \\
  -\cos\left(\phi + \psi \right) & -\sin\left( \phi + \psi \right) & 0 \\
\end{matrix}
\right)
\end{eqnarray*}
$$

\\
So for matrices near degenerate $\phi$ and $\psi$ can be expressed as:

$$ \begin{equation} \label{eq:msum}
\phi-\text{sgn}\left(m_{02}\right)\psi = \text{atan2}\left(m_{10},m_{20}\right)
\end{equation} $$

\\
and the equivalent quaternions respectively are:

$$
\frac{1}{\sqrt{2}}
  \left( \cos\left(\frac{\phi-\psi}{2}\right)  +
  \left( \sin\left(\frac{\phi-\psi}{2}\right),~
         \cos\left(\frac{\phi-\psi}{2}\right),~
        -\sin\left(\frac{\phi-\psi}{2}\right)\right) \right) \\
\frac{1}{\sqrt{2}}
  \left( \cos\left(\frac{\phi+\psi}{2}\right)  +
  \left( \sin\left(\frac{\phi+\psi}{2}\right),~
        -\cos\left(\frac{\phi+\psi}{2}\right),~
         \sin\left(\frac{\phi+\psi}{2}\right)\right) \right)
$$

\\
and $\phi$ and $\psi$ can be expressed as:

$$ \begin{equation} \label{qsum}
\phi+\text{sgn}\left(xz-wy\right)\psi = 2~\text{atan2}\left(x,w\right)
\end{equation} $$

\\
So my original implementation looked like the following and other versions I could found were pretty much equivalent.  Compute yaw by $\eqref{psiq2}$ and pitch by $\eqref{thetaq}$. Choose a degenerate zone threshold and use either $\eqref{phiq1}$ or $\eqref{qsum}$ for roll:


{% highlight c %}
void original(vec3_t* v, quat_t* q)
{
  float x=q->x, y=q->y, z=q->z, w=q->w;
  
  // half z-component of x' (negated)
  float xz = w*y-x*z;

  v->z = atan2f(x*y+w*z, .5f-(y*y+z*z));    // yaw
  v->y = atanf(xz/sqrtf(0.25f-xz*xz));      // pitch

  // roll computation broken in normal vs.
  // near degenerate case
  if (fabsf(xz) < YPR_GIMBAL_LOCK)
    v->x = atan2f(y*z+w*x, .5f-(x*x+y*y));
  else
    v->x = 2.f*atan2f(x,w) + sgn(xz)*v->z;
}
{% endhighlight %}


<br>

------

The duct-tape and super-glue
------

\\
Plugging the above into my quantization framework I noted something was rotten in Denmark in the form of unexpectly high errors.  Ripping out the quantization confirmed (peak error hitting ~5 degrees).  Changing the error plotting to $\abs{z}$ of $x'$ pointed the finger at the approaching denegerate case (suprise, suprise).  Sadly the error plot also showed that it was above my target peak error on over 15% of the range so the problem was beyond some minor tweak.

Fine.  First the denegerate zone test is based on an expression which is approaching one which can be changed to testing an expression approaching zero.  We also have some rounding problems such as subtract a small value from one.

For matrices we could compute $\theta$ from:

$$
\cos\left(\theta\right)^2 = m_{00}^2+m_{01}^2 = m_{12}^2+m_{22}^2
$$

and for quaternions the middle term expands to $\eqref{cos2_1}$, the right term to $\eqref{cos2_2}$ and they both reduce to $\eqref{cos2_3}$:

$$
\begin{eqnarray}
  \cos\left(\theta\right)^2 &=&  \left(\left(w^2-z^2\right)+\left(x^2-y^2\right)\right)^2 + ~\left(2xy+2wz\right)^2 \label{cos2_1} \\
                            &=&  ~(2wx+2yz)^2 + \left(\left(w^2-z^2\right)-\left(x^2-y^2\right)\right)^2  \label{cos2_2} \\
                            &=&  ((w+y)^2+(x-z)^2)((w-y)^2+(x+z)^2) \label{cos2_3}
							
\end{eqnarray}
$$

The infinite solution case is when $\cos\left(\theta\right)$ is zero and approaching zero can be used to detect the degenerate zone.  To state the obvious we can compute $\theta$ using the single parameter arc-tangent function:

$$
\theta = \text{atan}\left(\frac{\sin\left(\theta\right)}{\sqrt{\cos\left(\theta\right)^2}}\right)
$$


I'm targeting low hanging fruit, using $\eqref{cos2_1}$, carrying through some scale factors, pulling out some common expression and promoting to double computation yields:

{% highlight c %}
void revision_1(vec3_t* v, quat_t* q)
{
  double x=q->x, y=q->y, z=q->z, w=q->w;

  double t0 = (x+z)*(x-z);        // x^2-z^2
  double t1 = (w+y)*(w-y);        // w^2-y^2
  double xx = 0.5*(t0+t1);        // 1/2 x of x'
  double xy = x*y+w*z;            // 1/2 y of x'
  double xz = w*y-x*z;            // 1/2 z of x'
  double t  = xx*xx+xy*xy;        // cos(theta)^2
  double yz = 2.0*(y*z+w*x);      // z of y'

  v->z = (float)atan2(xy, xx);    // yaw   (psi)
  v->y = (float)atan(xz/sqrt(t)); // pitch (theta)

  if (t != 0) {
    v->x = (float)atan2(yz, t1-t0);
  else
    v->x = (float)(2.0*atan2(x,w) - sgnd(xz)*v->z);
}
{% endhighlight %}

Bingo! Well beyond my target.  Okay here are the error plots of some examples in the test framework:

* The original version: `orig`
* Simply promoting the original to double is useless: `orig_d`. Slightly better where is doesn't matter and the same error on the problem range.
* Change to doubles and new method: `rev_1`.  Not great accuracy but way beyond my target. Converting back to singles still has exploding error when $\abs{z}$ is approaching one. (not shown)
* Changing to singles for the majority of the range and using doubles for the remainder:  `rev_2 (5k)` is tight to my max peak error and `rev_2 (250k)` is roughly where the error start to otherwise explode (probability of ~.97 that the computation is in the common case of using singles).

<div id="error" style="width:100%"></div>

\\
<small>(EDIT: 20210529)</small> The original version of this post was mean and didn't have a work version in post (required you dig through the toy code and figure it out). Here is the *fma* based solution including required helper functions:

{% highlight c %}

inline float  sgn(float x)   { return copysignf(1.f,x); }

// ab+cd
inline float f32_mma(float a, float b, float c, float d)
{
  float t = c*d;
  float e = fmaf(c,d,-t);
  float f = fmaf(a,b, t);
  return f+e;
}

// ab-cd
inline float f32_mms(float a, float b, float c, float d)
{
  float t = c*d;
  float e = fmaf(c,d,-t);
  float f = fmaf(a,b,-t);
  return f-e;
}

// convert Tait-Bryan XYZ (commonly called Euler angles) to 
// a quaternion. 'revision_1_fma_s' in toy code
void xyz_to_quat(vec3_t* v, quat_t* q)
{
  float x=q->x, y=q->y, z=q->z, w=q->w;
  
  float t0 = (x+z)*(x-z);       // x^2-z^2
  float t1 = (w+y)*(w-y);       // w^2-y^2
  float xx = 0.5f*(t0+t1);
  float xy = f32_mma(x,y,w,z);
  float xz = f32_mms(w,y,x,z);
  float yz = 2.f*(f32_mma(y,z,w,x));
  float t  = xx*xx+xy*xy;

  v->z = atan2f(xy, xx);
  v->y = atanf(xz/sqrtf(t));

  if (t != 0)
    v->x = atan2f(yz, t1-t0);
  else
    v->x = 2.f*atan2f(x,w) - sgn(xz)*v->z;
}
{% endhighlight %}



\\
Extra credit:  the reverse transform is also a source of error.  My version is a direct translation of the math into code and isn't even bothering to pull out common sub-expressions that non-fast math optims cannot remove:

{% highlight c %}
void zyx_to_quat(quat_t* q, vec3_t* v)
{
  float hy = 0.5f*v->z;
  float hp = 0.5f*v->y;
  float hr = 0.5f*v->x;
  float ys = sinf(hy), yc = cosf(hy);
  float ps = sinf(hp), pc = cosf(hp);
  float rs = sinf(hr), rc = cosf(hr);

  q->w = rc*pc*yc + rs*ps*ys;
  q->x = rs*pc*yc - rc*ps*ys;
  q->y = rc*ps*yc + rs*pc*ys;
  q->z = rc*pc*ys - rs*ps*yc;
}
{% endhighlight %}

if this is modified to promote to double computation then the error drops a reasonable amount.  Example Using `rev_1` as the forward xform and double ZYX to quaternion gives us the error plot `rev1_id`.


<br>

------

Euler angles: ZXZ <small></small>
------

\\
Let's quickly walk through a proper Euler angle parameterization: ZXZ.  The quaternion rotation expression:

$$ 
\begin{eqnarray}
Q    &=& e^{\frac{\gamma}{2}\mathbf{z}}e^{\frac{\beta}{2}\mathbf{x}}e^{\frac{\alpha}{2}\mathbf{z}} \nonumber \\
     &=& \cos\left(\frac{\gamma+\alpha}{2}\right)\cos\left(\frac{\beta}{2}\right) + 
	 \left(
	  \cos\left(\frac{\gamma-\alpha}{2}\right)\sin\left(\frac{\beta}{2}\right),
	 ~\sin\left(\frac{\gamma-\alpha}{2}\right)\sin\left(\frac{\beta}{2}\right),
	 ~\sin\left(\frac{\gamma+\alpha}{2}\right)\cos\left(\frac{\beta}{2}\right)
	 \right) \label{zxz}
\end{eqnarray}
$$

\\
and using the similarity transform three times (tap, tap, tap) and converting into a matrix:

$$ 
\left(
\begin{array}{ccc}
 \cos (\alpha ) \cos (\gamma )-\sin (\alpha ) \cos (\beta ) \sin (\gamma ) & \sin (\alpha ) \cos (\beta )
   \cos (\gamma )+\cos (\alpha ) \sin (\gamma ) & \sin (\alpha ) \sin (\beta ) \\
 -\sin (\alpha ) \cos (\gamma )-\cos (\alpha ) \cos (\beta ) \sin (\gamma ) & \cos (\alpha ) \cos
   (\beta ) \cos (\gamma )-\sin (\alpha ) \sin (\gamma ) & \cos (\alpha ) \sin (\beta ) \\
 \sin (\beta ) \sin (\gamma ) & -\sin (\beta ) \cos (\gamma ) & \cos (\beta ) \\
\end{array}
\right)
$$

\\
For quaternions we don't need the matrix, directly from $\eqref{zxz}$: 

$$
\begin{eqnarray*}
 s &=& \frac{1}{2}\left(\gamma+\alpha\right) = \text{atan}\left(\frac{z}{w}\right) \\
 d &=& \frac{1}{2}\left(\gamma-\alpha\right) = \text{atan}\left(\frac{y}{x}\right) \\
 \\
 \gamma &=& s+d \\
 \alpha &=& s-d \\
 \\
 \beta &=& 2~\text{atan}\left(\sqrt{\frac{x^2+y^2}{w^2+z^2}}\right)
\end{eqnarray*}
$$

\\
Shoving these into a promote to double implementation gives:

{% highlight c %}
void quat_to_zxz(vec3_t* v, quat_t* q)
{
  double x=q->x, y=q->y, z=q->z, w=q->w;

  double s2 = x*x+y*y;      // sin(beta)^2
  double c2 = w*w+z*z;      // cos(beta)^2
  double s  = atan(z/w);    // (gamma+alpha)/2
  double d  = atan2(y,x);   // (gamma-alpha)/2
  v->x = s-d;               // alpha
  v->z = s+d;               // gamma

  if (c2 != 0.0)
    v->y = 2.0*atan(sqrt(s2/c2));
  else
    v->y = (0.5 > s2) ? 0 : PI;
}

{% endhighlight %}

\\
and a quick plot of above with single & double ZXZ to quaternion xform:

<div id="error2" style="width:100%"></div>

<script>

var xpos = new Array(100);

for (var i=0; i<100; i++) { xpos[i] = i*(1.0/99.0); };

var org = [0.000046, 0.000048, 0.000047, 0.000046, 0.000047, 0.000047, 0.000045, 0.000046, 0.000047, 0.000050, 0.000047, 0.000050, 0.000047, 0.000047, 0.000048, 0.000049, 0.000047, 0.000051, 0.000049, 0.000056, 0.000053, 0.000051, 0.000052, 0.000054, 0.000056, 0.000052, 0.000051, 0.000052, 0.000052, 0.000057, 0.000055, 0.000056, 0.000055, 0.000056, 0.000054, 0.000056, 0.000056, 0.000057, 0.000057, 0.000061, 0.000061, 0.000056, 0.000060, 0.000057, 0.000057, 0.000060, 0.000061, 0.000058, 0.000062, 0.000059, 0.000062, 0.000063, 0.000061, 0.000067, 0.000065, 0.000064, 0.000068, 0.000070, 0.000067, 0.000070, 0.000069, 0.000073, 0.000075, 0.000075, 0.000080, 0.000073, 0.000080, 0.000079, 0.000081, 0.000080, 0.000082, 0.000087, 0.000080, 0.000082, 0.000087, 0.000088, 0.000090, 0.000089, 0.000088, 0.000087, 0.000091, 0.000099, 0.000097, 0.000097, 0.000098, 0.000096, 0.000101, 0.000115, 0.000113, 0.000113, 0.000135, 0.000137, 0.000148, 0.000149, 0.000179, 0.000185, 0.000235, 0.000257, 0.000356, 4.945902];

var orgd = [0.000041, 0.000042, 0.000041, 0.000041, 0.000044, 0.000042, 0.000045, 0.000042, 0.000044, 0.000043, 0.000044, 0.000047, 0.000043, 0.000047, 0.000045, 0.000047, 0.000046, 0.000044, 0.000046, 0.000046, 0.000050, 0.000047, 0.000049, 0.000044, 0.000044, 0.000050, 0.000047, 0.000051, 0.000047, 0.000050, 0.000048, 0.000052, 0.000051, 0.000050, 0.000050, 0.000052, 0.000048, 0.000052, 0.000053, 0.000048, 0.000053, 0.000051, 0.000050, 0.000052, 0.000050, 0.000053, 0.000054, 0.000059, 0.000055, 0.000057, 0.000056, 0.000055, 0.000055, 0.000061, 0.000061, 0.000060, 0.000058, 0.000061, 0.000059, 0.000063, 0.000063, 0.000063, 0.000066, 0.000066, 0.000074, 0.000065, 0.000069, 0.000070, 0.000070, 0.000074, 0.000074, 0.000075, 0.000071, 0.000079, 0.000077, 0.000079, 0.000078, 0.000081, 0.000082, 0.000081, 0.000084, 0.000092, 0.000093, 0.000100, 0.000096, 0.000100, 0.000101, 0.000103, 0.000113, 0.000113, 0.000119, 0.000129, 0.000143, 0.000141, 0.000169, 0.000185, 0.000224, 0.000257, 0.000367, 4.9623];

var rev1 = [0.000022, 0.000022, 0.000022, 0.000022, 0.000022, 0.000022, 0.000022, 0.000022, 0.000024, 0.000023, 0.000023, 0.000023, 0.000022, 0.000022, 0.000023, 0.000023, 0.000022, 0.000024, 0.000022, 0.000022, 0.000022, 0.000023, 0.000023, 0.000022, 0.000025, 0.000022, 0.000021, 0.000021, 0.000021, 0.000021, 0.000023, 0.000021, 0.000021, 0.000023, 0.000023, 0.000022, 0.000022, 0.000021, 0.000023, 0.000022, 0.000024, 0.000023, 0.000024, 0.000023, 0.000023, 0.000023, 0.000024, 0.000024, 0.000023, 0.000024, 0.000024, 0.000024, 0.000023, 0.000024, 0.000024, 0.000024, 0.000025, 0.000024, 0.000025, 0.000023, 0.000029, 0.000028, 0.000028, 0.000029, 0.000029, 0.000029, 0.000029, 0.000027, 0.000029, 0.000028, 0.000028, 0.000029, 0.000029, 0.000028, 0.000030, 0.000030, 0.000028, 0.000026, 0.000026, 0.000027, 0.000027, 0.000028, 0.000027, 0.000027, 0.000028, 0.000027, 0.000027, 0.000027, 0.000028, 0.000028, 0.000027, 0.000029, 0.000029, 0.000028, 0.000028, 0.000029, 0.000031, 0.000029, 0.000030, 0.000031];

var rev1fma = [0.000019, 0.000020, 0.000019, 0.000020, 0.000019, 0.000019, 0.000020, 0.000019, 0.000020, 0.000020, 0.000020, 0.000020, 0.000020, 0.000019, 0.000020, 0.000019, 0.000019, 0.000020, 0.000020, 0.000020, 0.000021, 0.000020, 0.000020, 0.000020, 0.000020, 0.000020, 0.000020, 0.000021, 0.000023, 0.000020, 0.000022, 0.000020, 0.000021, 0.000020, 0.000020, 0.000021, 0.000020, 0.000021, 0.000020, 0.000021, 0.000021, 0.000020, 0.000021, 0.000022, 0.000020, 0.000020, 0.000021, 0.000022, 0.000021, 0.000021, 0.000022, 0.000022, 0.000023, 0.000022, 0.000023, 0.000023, 0.000022, 0.000021, 0.000022, 0.000022, 0.000026, 0.000027, 0.000028, 0.000027, 0.000027, 0.000027, 0.000028, 0.000027, 0.000028, 0.000028, 0.000027, 0.000027, 0.000029, 0.000027, 0.000027, 0.000028, 0.000027, 0.000026, 0.000026, 0.000026, 0.000026, 0.000026, 0.000026, 0.000026, 0.000027, 0.000025, 0.000026, 0.000026, 0.000026, 0.000026, 0.000027, 0.000027, 0.000027, 0.000028, 0.000027, 0.000026, 0.000027, 0.000027, 0.000029, 0.000028];

var rev2s5k =[0.000025, 0.000025, 0.000025, 0.000025, 0.000025, 0.000025, 0.000026, 0.000025, 0.000028, 0.000026, 0.000026, 0.000025, 0.000025, 0.000026, 0.000025, 0.000025, 0.000025, 0.000027, 0.000026, 0.000025, 0.000026, 0.000025, 0.000026, 0.000028, 0.000026, 0.000026, 0.000026, 0.000027, 0.000026, 0.000026, 0.000027, 0.000026, 0.000028, 0.000028, 0.000028, 0.000028, 0.000027, 0.000027, 0.000028, 0.000029, 0.000028, 0.000029, 0.000028, 0.000028, 0.000028, 0.000028, 0.000028, 0.000029, 0.000031, 0.000028, 0.000028, 0.000030, 0.000030, 0.000029, 0.000029, 0.000030, 0.000028, 0.000031, 0.000030, 0.000030, 0.000033, 0.000039, 0.000039, 0.000036, 0.000037, 0.000038, 0.000036, 0.000036, 0.000036, 0.000038, 0.000036, 0.000037, 0.000037, 0.000037, 0.000039, 0.000037, 0.000037, 0.000035, 0.000037, 0.000037, 0.000037, 0.000037, 0.000035, 0.000037, 0.000037, 0.000036, 0.000038, 0.000037, 0.000037, 0.000040, 0.000037, 0.000038, 0.000039, 0.000039, 0.000039, 0.000040, 0.000041, 0.000044, 0.000049, 0.000121];

var rev2s250k = [0.000025, 0.000025, 0.000025, 0.000025, 0.000025, 0.000025, 0.000026, 0.000025, 0.000028, 0.000026, 0.000026, 0.000025, 0.000025, 0.000026, 0.000025, 0.000025, 0.000025, 0.000027, 0.000026, 0.000025, 0.000026, 0.000025, 0.000026, 0.000028, 0.000026, 0.000026, 0.000026, 0.000027, 0.000026, 0.000026, 0.000027, 0.000026, 0.000028, 0.000028, 0.000028, 0.000028, 0.000027, 0.000027, 0.000028, 0.000029, 0.000028, 0.000029, 0.000028, 0.000028, 0.000028, 0.000028, 0.000028, 0.000029, 0.000031, 0.000028, 0.000028, 0.000030, 0.000030, 0.000029, 0.000029, 0.000030, 0.000028, 0.000031, 0.000030, 0.000030, 0.000033, 0.000039, 0.000039, 0.000036, 0.000037, 0.000038, 0.000036, 0.000036, 0.000036, 0.000038, 0.000036, 0.000037, 0.000037, 0.000037, 0.000039, 0.000037, 0.000037, 0.000035, 0.000037, 0.000037, 0.000037, 0.000037, 0.000035, 0.000037, 0.000037, 0.000036, 0.000038, 0.000037, 0.000037, 0.000040, 0.000037, 0.000038, 0.000039, 0.000039, 0.000039, 0.000040, 0.000041, 0.000029, 0.000030, 0.000031];

var rev1rd = [0.000012, 0.000013, 0.000013, 0.000012, 0.000012, 0.000013, 0.000013, 0.000013, 0.000014, 0.000013, 0.000012, 0.000013, 0.000013, 0.000014, 0.000014, 0.000013, 0.000013, 0.000013, 0.000014, 0.000014, 0.000014, 0.000013, 0.000014, 0.000014, 0.000013, 0.000013, 0.000014, 0.000014, 0.000014, 0.000014, 0.000014, 0.000014, 0.000014, 0.000014, 0.000014, 0.000014, 0.000014, 0.000015, 0.000014, 0.000014, 0.000014, 0.000014, 0.000014, 0.000015, 0.000015, 0.000014, 0.000015, 0.000014, 0.000015, 0.000015, 0.000015, 0.000014, 0.000015, 0.000014, 0.000015, 0.000015, 0.000014, 0.000014, 0.000014, 0.000015, 0.000015, 0.000015, 0.000015, 0.000015, 0.000015, 0.000015, 0.000015, 0.000015, 0.000015, 0.000016, 0.000015, 0.000016, 0.000016, 0.000015, 0.000015, 0.000015, 0.000015, 0.000018, 0.000018, 0.000018, 0.000018, 0.000018, 0.000018, 0.000018, 0.000018, 0.000018, 0.000018, 0.000018, 0.000018, 0.000018, 0.000019, 0.000018, 0.000018, 0.000018, 0.000018, 0.000019, 0.000019, 0.000019, 0.000019, 0.000019];

var euler1 = [0.000034, 0.000030, 0.000030, 0.000031, 0.000031, 0.000030, 0.000034, 0.000034, 0.000031, 0.000030, 0.000035, 0.000034, 0.000031, 0.000034, 0.000034, 0.000034, 0.000031, 0.000030, 0.000034, 0.000031, 0.000031, 0.000034, 0.000034, 0.000034, 0.000035, 0.000031, 0.000034, 0.000030, 0.000034, 0.000031, 0.000031, 0.000031, 0.000034, 0.000030, 0.000031, 0.000031, 0.000031, 0.000031, 0.000031, 0.000034, 0.000031, 0.000031, 0.000031, 0.000031, 0.000031, 0.000031, 0.000030, 0.000031, 0.000031, 0.000031, 0.000031, 0.000031, 0.000034, 0.000031, 0.000032, 0.000031, 0.000031, 0.000031, 0.000031, 0.000031, 0.000031, 0.000033, 0.000035, 0.000030, 0.000031, 0.000031, 0.000031, 0.000031, 0.000031, 0.000031, 0.000031, 0.000031, 0.000033, 0.000031, 0.000031, 0.000035, 0.000035, 0.000031, 0.000031, 0.000031, 0.000031, 0.000031, 0.000034, 0.000031, 0.000035, 0.000034, 0.000031, 0.000034, 0.000031, 0.000034, 0.000033, 0.000032, 0.000032, 0.000031, 0.000031, 0.000034, 0.000031, 0.000035, 0.000034, 0.000036];

var euler2 = [0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000017, 0.000019, 0.000017, 0.000017, 0.000019, 0.000019, 0.000019, 0.000020, 0.000019, 0.000020, 0.000019, 0.000019, 0.000020, 0.000019, 0.000019, 0.000020, 0.000020, 0.000020, 0.000020, 0.000020, 0.000020, 0.000020, 0.000020, 0.000020, 0.000020, 0.000020, 0.000021, 0.000021, 0.000020, 0.000021, 0.000021, 0.000021, 0.000021, 0.000021, 0.000020, 0.000021, 0.000021, 0.000021, 0.000021, 0.000021, 0.000021, 0.000021, 0.000021, 0.000021, 0.000021, 0.000021, 0.000021, 0.000021, 0.000023, 0.000021, 0.000021, 0.000021, 0.000021, 0.000021, 0.000021, 0.000021, 0.000021, 0.000021, 0.000023, 0.000021, 0.000021, 0.000022, 0.000021, 0.000022, 0.000021, 0.000021, 0.000021, 0.000022];

var etrace0 = { x: xpos, y: org,     name: 'orig'};
var etrace1 = { x: xpos, y: orgd,    name: 'orig_d'};
var etrace2 = { x: xpos, y: rev1,    name: 'rev_1'};
var etrace3 = { x: xpos, y: rev2s5k, name: 'rev_2s (5k)'};
var etrace4 = { x: xpos, y: rev2s250k, name: 'rev_2s (250k)'};
var etrace5 = { x: xpos, y: rev1rd,  name: 'rev_1id'};
var etrace6 = { x: xpos, y: rev1fma, name: 'rev_1_fma'};

var edata   = [etrace0,etrace1,etrace2,etrace3,etrace4,etrace5,etrace6];

var elayout = {
  xaxis: { title: 'abs(z)'},
  yaxis: { title: 'error', range:[0, 0.0001], hoverformat: 'g', exponentformat: 'power' }
};

Plotly.newPlot('error', edata, elayout);

var eulert1 = { x: xpos, y: euler1,  name: 'zxz'};
var eulert2 = { x: xpos, y: euler2,  name: 'zxz d'};
var edata2  = [eulert1,eulert2,etrace2,etrace5];
var elayout2 = {
  yaxis: { title: 'error', range:[0, 0.00005], hoverformat: 'g', exponentformat: 'power' }
};

Plotly.newPlot('error2', edata2, elayout2);


</script>
