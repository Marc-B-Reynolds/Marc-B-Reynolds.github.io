---
layout:       post
title:        On quaternion/rotation matrix conversions and errors
categories:   [quaternions]
tags:         [rotations, conversion]
plotly:       true
disqus:       false
description:  Skims the math of quaternions to rotations matrices and back with empirical errors including small angle special case.
---

Toy code can be found: [HERE](http://github.com/Marc-B-Reynolds/Stand-alone-junk/blob/master/src/Posts/q2mat.c).  All reported values are using clang without enabling floating-point contractions.
<br>

------

Prelim <small>the dead horse beating</small>
------

\\
Skimming the math.  Starting with a quaternion $Q=w+\left(x,y,z\right)$ then we can rotate $\mathbf{p}$ by:

$$ \begin{equation} \label{eq:sxform}
\mathbf{p}' = Q\mathbf{p}Q^{-1}
\end{equation} $$

\\
and if $Q$ is unit magnitude this reduces to:

$$ \begin{equation} \label{eq:rsxform}
\mathbf{p}' = Q\mathbf{p}Q^*
\end{equation} $$

\\
If $Q$ with magitude of $s$ is transformed by $\eqref{eq:rsxform}$ then the result is the composition of the rotation and uniform scaling of $s^2$.

To create a matrix we need to apply the rotation to the basis set to form our three equations:

$$ 
\begin{eqnarray*}
\mathbf{x} = \left(1,0,0\right) \\
\mathbf{y} = \left(0,1,0\right) \\
\mathbf{z} = \left(0,0,1\right) \\
\end{eqnarray*}
$$

which expanded and reduced gives:

$$ \begin{eqnarray}
\mathbf{x}' & = & Q\mathbf{x}Q^{*} \nonumber \\
            & = & \left(w^2 + x^2 - y^2 - z^2, ~2\left(xy+wz\right), ~2\left(xz-wy\right) \right) \label{xp}    \\
            & = & \left(1 - 2 \left(y^2+z^2\right), ~2\left(xy+wz\right), ~2\left(xz-wy\right) \right) \label{xpr} \\
			\nonumber \\
\mathbf{y}' & = & Q\mathbf{y}Q^{*} \nonumber \\
            & = & \left(2 \left(xy-wz\right), ~w^2 - x^2 + y^2 - z^2, ~2\left(wx+yz\right) \right)    \label{yp} \\
            & = & \left(2 \left(xy-wz\right), ~1 - 2\left(x^2+z^2\right), ~2\left(wx+yz\right) \right) \label{ypr} \\
			\nonumber \\
\mathbf{z}' & = & Q\mathbf{z}Q^{*} \nonumber \\
            & = & \left(2 \left(wy+xz\right), ~2\left(yz-wx\right), ~w^2 - x^2 - y^2 + z^2\right) \label{zp} \\
            & = & \left(2 \left(wy+xz\right), ~2\left(yz-wx\right), ~1 - 2\left(x^2+y^2\right) \right) \label{zpr} 
\end{eqnarray} $$

<br>

------

Quaternion to rotation matrix
------

\\
Sticking to the math convention of column vectors, then we can shove the (not reduced) equations $\eqref{xp} \eqref{yp} \eqref{zp}$ into the columns giving:

$$ \begin{equation} \label{2matndr}
{
\left( \begin{array}{ccc}
w^2 + x^2 - y^2 - z^2 & 2\left(xy-wz\right)       &  2\left(xz+wy\right)       \\
2 \left(xy+wz\right)       & w^2 - x^2 + y^2 - z^2 &  2\left(yz-wx\right)      \\
2 \left(xz-wy\right)       & 2\left(yz+wx\right)       & w^2 - x^2 - y^2 + z^2
\end{array} \right)
} \end{equation}
$$

\\
and we can reduce the diagonal elements since $Q\cdot Q = 1$ as per $\eqref{xpr} \eqref{ypr} \eqref{zpr}$:

$$ \begin{equation} \label{2matstd}
{
\left( \begin{array}{ccc}
1 - 2 \left(y^2+z^2\right) & 2\left(xy-wz\right)       &  2\left(xz+wy\right)   \\
2 \left(xy+wz\right)       & 1 - 2\left(x^2+z^2\right) &  2\left(yz-wx\right)   \\
2 \left(xz-wy\right)       & 2\left(yz+wx\right)       & 1 - 2\left(x^2+y^2\right) 
\end{array} \right)
} \end{equation}
$$

\\
  Translating both into C gives:

{% highlight c %}

// not reducing the diagonal. adds ~4 (basic simple scalar) ops
// vs. the standard version.
void quat_to_mat33_ndr(mat33_t* m, quat_t* q)
{
  float x  = q->x, y  = q->y, z  = q->z, w  = q->w;
  float xx = x*x,  yy = y*y,  zz = z*z,  ww = w*w;
  float tx = 2*x,  ty = 2*y,  tz = 2*z;
  float xy = ty*x, xz = tz*x, yz = ty*z;
  float wx = tx*w, wy = ty*w, wz = tz*w;
  float t0 = ww-zz;
  float t1 = xx-yy;
  
  m->m00 = t0+t1;
  m->m11 = t0-t1;
  m->m22 = ww-xx-yy+zz;

  m->m10 = xy+wz; m->m01 = xy-wz;
  m->m20 = xz-wy; m->m02 = xz+wy;
  m->m21 = yz+wx; m->m12 = yz-wx;
}

// commonly seen variant
void quat_to_mat33_std(mat33_t* m, quat_t* q)
{
  float x  = q->x, y  = q->y, z  = q->z, w  = q->w;
  float tx = 2 *x, ty = 2 *y, tz = 2 *z;
  float xx = tx*x, yy = ty*y, zz = tz*z;
  float xy = ty*x, xz = tz*x, yz = ty*z;
  float wx = tx*w, wy = ty*w, wz = tz*w;

  m->m00 = 1.f-(yy+zz); m->m11 = 1.f-(xx+zz); m->m22 = 1.f-(xx+yy);

  m->m10 = xy+wz; m->m01 = xy-wz;
  m->m20 = xz-wy; m->m02 = xz+wy;
  m->m21 = yz+wx; m->m12 = yz-wx;
}

{% endhighlight %}

<br>

<div class="alert alert-info" role="alert">
<strong>NOTE:</strong> An identical problem is extracting one or more axes from a quaterion.  Say $\mathbf{x}, \mathbf{y}, \mathbf{z}$ respectively represent: right, forward, up. We can find the transformed direction(s) using the corresponding primed equation(s).
<br>

<br>
The reason for this aside that because I've sometime seen code that extracts two of the directions as above and builds the third using a cross product.  For the cross product method to be a performance win requires that it ends up being faster than three additions.
</div>

Another aside: In cases where input quaternions are moving away from unity due to compounding errors it is possible to fold the corrective step into the transform instead of renormalizing as a pre-step.

<br>

------

Rotation matrix to quaternion
------

\\
Let's rename the elements of $\eqref{2matndr}$ as:

$$
{
\left( \begin{array}{ccc}
m_{00} & m_{01} & m_{02} \\
m_{10} & m_{11} & m_{12} \\
m_{20} & m_{21} & m_{22}
\end{array} \right)
}
$$

\\
then the sum and difference of off diagonal pairs give:

$$
\begin{eqnarray}
  m_{10} + m_{01} & = & 4xy \label{4xy} \\
  m_{10} - m_{01} & = & 4wz \label{4wz} \\
  m_{02} + m_{20} & = & 4xz \label{4xz} \\
  m_{02} - m_{20} & = & 4wy \label{4wy} \\
  m_{21} + m_{12} & = & 4yz \label{4yz} \\
  m_{21} - m_{12} & = & 4wx \label{4wx}
\end{eqnarray}
$$


\\
and the permutations of *trace* like functions:

$$
\begin{eqnarray*}
  m_{00} + m_{11} + m_{22} & = & 3 w^2 - \left(x^2 + y^2 + z^2 \right)  & = & -1+4w^2 \\ 
  m_{00} + m_{11} - m_{22} & = & \left(w^2 + x^2 + z^2 \right)- 3 y^2   & = & 1-4y^2  \\ 
  m_{00} - m_{11} + m_{22} & = & \left(w^2 + x^2 + y^2 \right) - 3 z^2  & = & 1-4z^2  \\ 
  m_{00} - m_{11} - m_{22} & = & 3 x^2 - \left(w^2 + y^2 + z^2\right)  & = & -1+4x^2 \\ 
\end{eqnarray*}
$$

\\
which can be reworked to:

$$
\begin{eqnarray}
  4w^2 & = & 1 + m_{00} + m_{11} + m_{22} \label{4w2} \\
  4x^2 & = & 1 + m_{00} - m_{11} - m_{22} \label{4x2}\\
  4y^2 & = & 1 - m_{00} + m_{11} - m_{22} \label{4y2}\\ 
  4z^2 & = & 1 - m_{00} - m_{11} + m_{22} \label{4z2}
\end{eqnarray}
$$

\\
If we know or determine an element has a *sufficiently* large magnitude we can choose the equation for that element from the second set of equations and three equations from the first set which are products of the that and the remaining elements.

As an example take a special case version where we know that the magnitude of the angle is not greater than $\frac{\pi}{2}$ which in turn ensures that the scalar $(w)$ is the largest element.  Then we need $\eqref{4w2}$ from the second set and $\eqref{4wx} \eqref{4wy} \eqref{4wz}$ from the first.

<br>

{% highlight c %}
void mat33_2_quat_small(quat_t* q, mat33_t* m)
{
  float m00=m->m00, m01=m->m01, m02=m->m02;
  float m10=m->m10, m11=m->m11, m12=m->m12;
  float m20=m->m20, m21=m->m21, m22=m->m22;
  float t  = 1 + m00 + m11 + m22;  // 4w^2
  float s  = 0.5f*rsqrtf(t);       // 1/(4w)
  quat_set(q, s*(m21-m12), s*(m02-m20), s*(m10-m01), s*t);
}
{% endhighlight %}

<br>

#### Standard construction

The standard construction for the general case is based on the observation that since we're producing a unit quaternion, then the smallest that the largest magnitude component can be is $\frac{1}{2}$. Given that we can simply walk through the components and choose the first that's at least that large.

{% highlight c %}
void mat33_to_quat_std(quat_t* q, mat33_t* m)
{
  float m00=m->m00, m01=m->m01, m02=m->m02;
  float m10=m->m10, m11=m->m11, m12=m->m12;
  float m20=m->m20, m21=m->m21, m22=m->m22;
  float t1 = m22+m11;
  float t0 = m00+t1;

  if (t0 > 0.f) {
    float w = 1.f+t0, s = 0.5f*rsqrtf(w);
    quat_set(q, s*(m21-m12), s*(m02-m20), s*(m10-m01), s*w); // w
    return;
  }

  t0 = m00-t1;
  
  if (t0 > 0.f) {
    float x = 1.f+t0, s = 0.5f*rsqrtf(x);
    quat_set(q, s*x, s*(m01+m10), s*(m02+m20), s*(m21-m12)); // x
    return;
  }

  t0 = m11-m22; t1 = 1.f-m00;
  
  if (t0 > 0.f) {
    float y = t1+t0, s = 0.5f*rsqrtf(y);
    quat_set(q, s*(m01+m10), s*y, s*(m12+m21), s*(m02-m20)); // y
  }
  else {
    float z = t1-t0, s = 0.5f*rsqrtf(z);
    quat_set(q, s*(m02+m20), s*(m12+m21), s*z, s*(m10-m01)); // z
  }
}
{% endhighlight %}

\\
Given uniform random input then the (in code order) probabilities are: $\frac{2}{3},\frac{1}{9},\frac{1}{9},\frac{1}{9} $  The first from $\frac{2}{\pi}\text{acos}\left(\frac{1}{2}\right) = \frac{2}{3}$.

<br>

#### Mike Day's construction

Insomniac Games (Mike Day) presented a paper[^day] in 2015 which makes the following observations. The sign of diagonal components imply $\eqref{2matstd}$:

$$ 
\begin{eqnarray*}
m_{00} < 0 \implies y^2+z^2 > \frac{1}{2} \\
m_{11} < 0 \implies x^2+z^2 > \frac{1}{2} \\
m_{22} < 0 \implies x^2+y^2 > \frac{1}{2}
\end{eqnarray*}
$$

\\
and the permutations of sum/difference of pairs of diagonal elements imply:

$$
\begin{eqnarray*}
m_{00}+m_{11} < 0 \implies w^2 < z^2 \\
m_{00}-m_{11} < 0 \implies x^2 < y^2 \\
m_{00}+m_{22} < 0 \implies w^2 < y^2 \\
m_{00}-m_{22} < 0 \implies x^2 < z^2 \\ 
m_{11}+m_{22} < 0 \implies w^2 < x^2 \\
m_{11}-m_{22} < 0 \implies y^2 < z^2
\end{eqnarray*}
$$

They use the above to form a method that insures to choose a sufficiently large element.  The following is a slightly massaged version of Day:

{% highlight c %}
void mat33_2_quat_day(quat_t* q, mat33_t* m)
{
  float m00=m->m00, m01=m->m01, m02=m->m02;
  float m10=m->m10, m11=m->m11, m12=m->m12;
  float m20=m->m20, m21=m->m21, m22=m->m22;
  float e0;

  if (m22 >= 0) {
    float a = m00+m11;
    float b = m10-m01;
    float c = 1.f+m22;

    if (a >= 0) { 
      e0 = c + a;
      quat_set(q, m21-m12, m02-m20, b, e0); // w
    }
    else { 
      e0 = c - a;
      quat_set(q, m02+m20, m21+m12, e0, b); // z
    }
  }
  else {
    float a = m00-m11;
    float b = m10+m01;
    float c = 1.f-m22;
    
    if (a >= 0) {
      e0 = c + a;
      quat_set(q, e0, b, m02+m20, m21-m12); // x
    }
    else {
      e0 = c - a;
      quat_set(q, b, e0, m21+m12, m02-m20); // y
    }
  }

  quat_scale(q, 0.5f*rsqrtf(e0));
}
{% endhighlight %}

\\
Given uniform random input then the outer and both inner probabilities are: $\frac{1}{2}$.

<br>

#### Angle errors

\\
Let's examine round-trip errors.  We generate a uniform quaternion $A$, convert to a matrix, back to a quaternion $B$ and then find the relative quaternion $X$:

$$ X = BA^* = \mathbf{v} + w $$

\\
then measure the represented angle, multiple by two (to get 3D rotation measure instead of in $\mathbb{H}$) and convert to degrees:

$$ \frac{360}{\pi}\text{atan}\left(\frac{\left\Vert \mathbf{v} \right\Vert }{\abs{w}}\right)$$

\\
This choice of measure ignores the magnitudes of $A$ and $B$. The trace names are to-matrix-method/to-quat-method.  The x-axis is normalized.

<p align="center"><div id="fig1" style="width:100%"></div></p>

\\
There's a jump in trace behavior at $\frac{2}{3}$ because our x-axis is a function of $w$ and this the point where $w=\frac{1}{2}$. The forward transform without diagonal reduction has lower peak error, but *potentially* more interesting is the fact that this method is better at angle-preserving (conformal) than the standard.


### Variants <small>you probably don't care</small>

<br>

#### Branch free Day
\\
The reworking of Day above is to point out the we can formulate a branch-free version. Let's call a function to isolate the sign bit of a float $sb$ and denote an integer XOR as $\oplus$. The initial branch is selected based on the sign of $m_{22}$. In both cases we're computing three temp variables $a,b,c$ which can be expressed in terms of that sign: $s_0 = \text{sb}\left(m_{22}\right)$.  Then inside the inner four cases we're computing $e_0$ (the chosen sufficiently large element) which can be expressed in terms of the sign of $a$: $s_1 = \text{sb}\left(a\right)$. ($c$ and $e_0$ can use abs instead).  The remaining two elements (unnamed temp expressions) can be formed using $s_1$ and the parity of two sign values:  $s_p = s_0 \oplus s_1$.  If we rename $b$ to $e_1$, then so far we have:

$$ \begin{eqnarray*}
s_0 & = & \text{sb}\left(m_{22}\right)  \\
a   & = & m_{00} + \left(m_{11} \oplus s_0\right) \\
c   & = & 1      + \left(m_{22} \oplus s_0\right) = 1 + \abs{m_{22}} \\
s_1 & = & \text{sb}\left(a\right) \\
s_p & = & s_0 \oplus s_1 \\
e_0 & = & c + \left(a \oplus s_1\right) = c + \abs{a} \\
e_1 & = & m_{10} - \left(m_{01} \oplus s_0\right)  \\
e_2 & = & m_{02} - \left(m_{20} \oplus s_p\right) \\
e_3 & = & m_{21} - \left(m_{12} \oplus s_1\right)
\end{eqnarray*} $$

\\
To complete the component computation we simply need to scale all four by $\frac{1}{2 \sqrt{e_0}}$.  This leave placing them in memory, if we has storage order of $\left(x,~y,~z,~w\right)$ and think of our elements as the ordered list: $\left(e_3,~e_2,~e_1,~e_0 \right) $, then we apply the follow two rules (input to the second is renaming the output of the first)

* if $s_0$ is set $\left(e_3,~e_2,~e_1,~e_0 \right) \rightarrow \left( e_1,~e_0,~e_3,~e_2 \right) $, so flip the first and last two. 
* if $s_p$ is set $\left(e_3,~e_2,~e_1,~e_0 \right) \rightarrow \left( e_2,~e_3,~e_0,~e_1 \right) $, so flip the first and second pairs. 

Alternately we can simply compute a value to XOR as below:

{% highlight c %}
void mat33_to_quat_day_bf(quat_t* q, mat33_t* m)
{
  float m00=m->m00, m01=m->m01, m02=m->m02;
  float m10=m->m10, m11=m->m11, m12=m->m12;
  float m20=m->m20, m21=m->m21, m22=m->m22;

  // isolate a sufficently large component
  uint32_t s0 = sb(m22);
  float    a  = m00 + fxor(m11,s0);              // m00 + (m11^s0)
  uint32_t s1 = sb(a);
  uint32_t sp = s0^s1;
  float    e0 = 1.f + fxor(m22,s0) + fxor(a,s1); // 1 + |m22| + |a|
  float    s  = 0.5f*rsqrtf(e0);

  // the other three components
  float    e1 = m10-fxor(m01,s0);                // m10 - (m01^s0)
  float    e2 = m02-fxor(m20,sp);                // m02 - (m02^sp)
  float    e3 = m21-fxor(m12,s1);                // m21 - (m12^s1)

  e0 *= s; e1 *= s; e2 *= s; e3 *= s;            // scaling

  uint32_t id = (s0 >> 30)^(sp >> 31);           // write ordering

  q->f[3^id] = e0;
  q->f[2^id] = e1;
  q->f[1^id] = e2;
  q->f[0^id] = e3;
}
{% endhighlight %}

<br>


#### Branch free <small>higher error</small>

\\
Another branch-free possibility is to use the second set of equations to compute the magnitudes of all the elements, choose the scalar to be positive and compute the sign of the bivector components by the first set of equations.

{% highlight c %}
void mat33_2_quat_bf0(quat_t* q, mat33_t* m)
{	
  float m00=m->m00, m01=m->m01, m02=m->m02;
  float m10=m->m10, m11=m->m11, m12=m->m12;
  float m20=m->m20, m21=m->m21, m22=m->m22;

  float a = 1.f+m00;
  float b = 1.f-m00;
  float c = m11+m22;
  float d = m11-m22;

  // abs to prevent errors, could use max against zero
  // or a bias amount as per below.
  float w = fabsf(a+c); // 4w^2
  float x = fabsf(a-c); // 4x^2
  float y = fabsf(b+d); // 4y^2
  float z = fabsf(b-d); // 4z^2

  // assumes rsqrtf(0) = inf, otherwise needs to add a bias
  w *= 0.5f*rsqrtf(w);
  x *= 0.5f*rsqrtf(x);
  y *= 0.5f*rsqrtf(y);
  z *= 0.5f*rsqrtf(z);

  // use first set for the signs
  x = copysignf(x, m21-m12);
  y = copysignf(y, m02-m20);
  z = copysignf(z, m10-m01);

  quat_set(q, x,y,z,w);
}
{% endhighlight %}

\\
which as promised results in much higher errors:

<div id="fig2" style="width:100%"></div>

\\
Although there are cases where these errors are low enough to be usable the performace is worse that some better methods...


<br>

#### Almost branch free

\\
Yet another option is to notice that special case *small angle* version performs rather well.  These plots again are:

<div id="fig3" style="width:100%"></div>

\\
So we can formulate a version that has a highly predictable branch that handles the common case and otherwise fall-back.

{% highlight c %}
void mat33_to_quat_bf1(quat_t* q, mat33_t* m)
{
  float m00=m->m00, m01=m->m01, m02=m->m02;
  float m10=m->m10, m11=m->m11, m12=m->m12;
  float m20=m->m20, m21=m->m21, m22=m->m22;
  float t  = 1.f + m00 + m11 + m22;   // 4w^2

  if (t > BF1_CUT) {
    float s  = 0.5f*rsqrtf(t);
    quat_set(q, s*(m21-m12), s*(m02-m20), s*(m10-m01), s*t);
    return;
  }

  // add handling other cases here
}

{% endhighlight %}

\\
A lazy table of empirical numbers:  some ad-hoc cut points, probability and round-trip errors using the standard method of quat-to-matrix:

{: .center }
| BF1_CUT |probability| rt-error |
| ------  | --------- | -------  |
| 0.00001 | 0.998002  | 0.042664 |
| 0.0001  | 0.993634  | 0.014877 |
| 0.0002  | 0.990994  | 0.010233 |
| 0.0004  | 0.987265  | 0.007769 |
| 0.0008  | 0.981993  | 0.004991 |
| 0.0016  | 0.974531  | 0.003587 |
| 0.0032  | 0.963969  | 0.002935 |
| 0.0064  | 0.949045  | 0.001946 |
| 0.0128  | 0.927966  | 0.001367 |
| 0.0256  | 0.898184  | 0.000975 |
| 0.0512  | 0.856216  | 0.000753 |
| 0.1024  | 0.797129  | 0.000528 |

<br>


------

Overkill <small>but lower error</small>
------
{:#lowerror}

\\
Let's run through a method which includes every component of the rotation matrix in each component of the resulting quaternion. The point of that exercise is a light filtering of noise present in the input.  First we'll rewrite the rotation matrix as the result of converting an axis-angle conversion with unit vector $(x,~y,~z)$ and angle $\theta \in [0,\pi]$:

$$
\left(
\begin{array}{ccc}
 \cos \theta +x^2 (1-\cos \theta ) & x y (1-\cos \theta )-z \sin \theta  & x z (1-\cos \theta)+y \sin \theta  \\
 x y (1-\cos \theta )+z \sin \theta  & \cos \theta +y^2 (1-\cos \theta ) & y z (1-\cos \theta )-x \sin \theta  \\
 x z (1-\cos \theta )-y \sin \theta  & y z (1-\cos \theta) + x \sin \theta & \cos \theta +z^2 (1-\cos \theta ) \\
\end{array}
\right)
$$

\\
**NOTE:** everywhere else I'm using $x,y,z$ for quaternion components and here for the unit vector of the axis angle.

Expanding and naming off diagonal differences:

$$
\begin{eqnarray*}
  s_1 = m_{21} - m_{12} = 2x \sin(\theta) \\
  s_2 = m_{02} - m_{20} = 2y \sin(\theta) \\
  s_3 = m_{01} - m_{10} = 2z \sin(\theta) 
\end{eqnarray*}
$$

\\
and sums:

$$
\begin{eqnarray*}
  a_1 = m_{21} + m_{12} = 2yz \left(1-\cos \theta \right) \\
  a_2 = m_{02} + m_{20} = 2xz \left(1-\cos \theta \right) \\
  a_3 = m_{01} + m_{10} = 2xy \left(1-\cos \theta \right) 
\end{eqnarray*}
$$

\\
and (what I'm calling) trace like:

$$
\begin{array}{lll}
  t_0 = & m_{00} + m_{11} + m_{22} = & 1 + 2~\cos \theta \\
  t_1 = & m_{00} - m_{11} - m_{22} = & 2~x^2 \left(1-\cos \theta \right)-1 \\
  t_2 = & m_{11} - m_{00} - m_{22} = & 2~y^2 \left(1-\cos \theta \right)-1 \\
  t_2 = & m_{22} - m_{00} - m_{11} = & 2~z^2 \left(1-\cos \theta \right)-1
\end{array}
$$

\\
We can combine these as follows:

$$
\begin{array}{lll}
  r_0 = (t_0+1)^2 + s_1^2 + s_2^2 + s_3^2 = & 8\left(1+\cos \theta \right) \\
  r_1 = (t_1+1)^2 + s_1^2 + a_2^2 + a_3^2 = & 8\left(x^2 \left(1-\cos \theta \right)\right) \\
  r_2 = (t_2+1)^2 + s_2^2 + a_1^2 + a_3^2 = & 8\left(y^2 \left(1-\cos \theta \right)\right) \\
  r_3 = (t_3+1)^2 + s_3^2 + a_1^2 + a_2^2 = & 8\left(z^2 \left(1-\cos \theta \right)\right) 
\end{array}
$$

\\
a simple rewrite to make the [half-angle](https://en.wikipedia.org/wiki/List_of_trigonometric_identities#Multiple-angle_formulae) identities jump out:

$$
\begin{array}{lll}
  r_0 = 16 \left(\frac{1+\cos \theta}{2} \right) \\
  r_1 = 16~x^2 \left(\frac{1-\cos \theta}{2} \right) \\
  r_2 = 16~y^2 \left(\frac{1-\cos \theta}{2} \right) \\
  r_3 = 16~z^2 \left(\frac{1-\cos \theta}{2} \right) 
\end{array}
$$

\\
So we can compute the magnitudes of the quaternion components as follows:

$$
\begin{array}{lll}
  |q_w| = \frac{\sqrt{r_0}}{4} \\
  |q_x| = \frac{\sqrt{r_1}}{4} \\
  |q_y| = \frac{\sqrt{r_2}}{4} \\
  |q_z| = \frac{\sqrt{r_3}}{4} 
\end{array}
$$

\\
The remaining issue is then determining the signs of the components and the $s_i$ terms have that information. Plot the error of that implementation (vs. Day) would give us:

<div id="fig4" style="width:100%"></div>

\\
So it all falls apart when approaching the maximal rotation angle. The problem is $\sin \theta$ is approaching zero. We can correct this by noting the most important information is in the largest magnitude component $q_i$ and we could choose that to be the one that's positive (or equivalently known sign). Specifically the $a_i$ term which includes the known and desired sign components. With this correction we have:

<div id="fig5" style="width:100%"></div>

\\
Let's do some points:

* *dp* is promoting input to doubles, computing and demoting to single. For single precision this is less work than being fancy. If working precision is doubles you'd have to go with extended precision methods (*TwoProduct*, compensated summations, etc.) to squeeze out the additonal error reduction.
* Although it's structure helps with noise filtering, performing the standard quaternion to matrix method (reducing the diagaonals) is still a major poke in the eye. Just hide all the *std* methods now..thanks!
* All of this added computational complexity does buy us anything for small rotation angles (vs. Day & "Almost branchfree"..provided NDR quat-to-mat)

\\
Since NDR is the most effective error reduction, let's do a spitball of it promoted to double computation. (reminder: the toy code will let you play with whatever combo you want)

<div id="fig7" style="width:100%"></div>


\\
Now we're not really attempting to filter a noisy input we can relax to ensuring $w$ is sufficiently large to not increase the error in the case we're interested it to not have a branch misprediction nightmare. (Going branchfree is possible and would be the "thing" todo in a SIMD implementation. For scalar it doesn't seem worthwhile). Making the simplifing assumption that rotations are uniformly random distributed then the normalized [cumulative volume](http://marc-b-reynolds.github.io/quaternions/2017/11/10/AveRandomRot.html) (cumulative/total) of a rotation is:

$$ \frac{\theta - \sin \theta}{\pi} $$

\\
for a 3D rotation angle of $\theta$. Plugging in how the quaternion scalar $w$ is related to $\theta = 2 ~\text{acos} (w)$ yields:

$$ \frac{2 ~\text{acos} (w) - 2w\sqrt{1-w^2}}{\pi} $$

which is pretty close to a linear decreasing function:

<div id="fig6" style="width:100%"></div>

\\
So if the first test is `w > CUT` then the probablity of that block is $w$ fed into the above expression (assuming uniform input and better if angles are skewed toward small).  An implementation:


{% highlight c %}
// example cut points:
// @ 0.1015625 reaches 'w' case 87.1% <- this is playin' it safe 
// @ 0.0678711 reaches 'w' case 91.3% <- this is close to min
//					 value w/o error increase

// example choice
#define OKC_CUT 0.0678711f

void mat33_to_quat_okc(quat_t* q, mat33_t* m)
{	
  float m00=m->m00, m01=m->m01, m02=m->m02;
  float m10=m->m10, m11=m->m11, m12=m->m12;
  float m20=m->m20, m21=m->m21, m22=m->m22;

  float t0 = m00+m11+m22, d0 = t0+1.f;
  float t1 = m00-m11-m22, d1 = t1+1.f;
  float t2 = m11-m00-m22, d2 = t2+1.f;
  float t3 = m22-m00-m11, d3 = t3+1.f;
  
  float s1 = m21-m12, s2 = m02-m20, s3 = m10-m01;
  float a1 = m21+m12, a2 = m02+m20, a3 = m10+m01;
  
  float w = 0.25f*sqrtf(d0*d0 + s1*s1 + s2*s2 + s3*s3);
  float x = 0.25f*sqrtf(d1*d1 + s1*s1 + a2*a2 + a3*a3);
  float y = 0.25f*sqrtf(d2*d2 + s2*s2 + a1*a1 + a3*a3);
  float z = 0.25f*sqrtf(d3*d3 + s3*s3 + a1*a1 + a2*a2);

  // could be make branchfree (thinkin' SIMD here)
  if (w > OKC_CUT) {
    x = copysignf(x, s1);
    y = copysignf(y, s2);
    z = copysignf(z, s3);
  }
  else {
    if (x > 0.25f) {
      y = copysignf(y,a3);
      z = copysignf(z,a2);
      w = copysignf(w,s1);
    }
    else if (y > 0.25f) {
      x = copysignf(x, a3);
      z = copysignf(z, a1);
      w = copysignf(w, s2);
    }
    else {
      x = copysignf(x, a2);
      y = copysignf(y, a1);
      w = copysignf(w, s3);
    }
  }
  quat_set(q, (float)x,(float)y,(float)z,(float)w);
}
{% endhighlight %}

<br>


------

References and Footnotes
------

[^day]: **"Converting a Rotation Matrix to a Quaternion"**, Mike Day, 2015 ([download page](http://www.insomniacgames.com/converting-a-rotation-matrix-to-a-quaternion/))

<script>

var bax = new Array(80);

{ for(var i=0; i<80; i++) bax[i] = i/79.0; }

const stdstd = [3.366163e-07, 7.089606e-07, 1.002567e-06, 1.467112e-06, 1.965123e-06, 2.060151e-06, 2.680844e-06, 3.431232e-06, 3.765651e-06, 3.850529e-06, 4.297061e-06, 5.300426e-06, 5.037166e-06, 6.900605e-06, 7.025214e-06, 8.815179e-06, 8.059801e-06, 9.133669e-06, 1.013469e-05, 1.025556e-05, 1.104310e-05, 1.280408e-05, 1.284779e-05, 1.341579e-05, 1.380421e-05, 1.606909e-05, 1.868217e-05, 2.235794e-05, 2.254608e-05, 2.664295e-05, 2.799844e-05, 2.702300e-05, 2.964451e-05, 3.336401e-05, 3.319995e-05, 3.651717e-05, 3.869078e-05, 3.879974e-05, 4.364011e-05, 4.584912e-05, 5.431133e-05, 6.020936e-05, 6.017704e-05, 6.368387e-05, 6.587051e-05, 7.539200e-05, 8.067569e-05, 8.180513e-05, 8.839987e-05, 9.791440e-05, 1.004588e-04, 1.025081e-04, 1.215223e-04, 1.095087e-04, 1.541106e-05, 1.501028e-05, 1.511791e-05, 1.578418e-05, 1.541071e-05, 1.514047e-05, 1.575093e-05, 1.585129e-05, 1.577171e-05, 1.574984e-05, 1.596803e-05, 1.575143e-05, 1.604615e-05, 1.594718e-05, 1.552819e-05, 1.609223e-05, 1.552909e-05, 1.577294e-05, 1.581677e-05, 1.519899e-05, 1.591611e-05, 1.574458e-05, 1.565862e-05, 1.579813e-05, 1.572727e-05, 1.575847e-05];

const ndrstd = [2.444860e-07, 5.222806e-07, 7.442844e-07, 1.006310e-06, 1.459467e-06, 1.721504e-06, 2.588514e-06, 2.920642e-06, 3.033302e-06, 3.556354e-06, 3.844039e-06, 4.009854e-06, 4.704491e-06, 5.158537e-06, 5.627595e-06, 6.705708e-06, 7.230056e-06, 7.696886e-06, 7.819148e-06, 7.650115e-06, 8.471955e-06, 8.790127e-06, 9.268349e-06, 1.038301e-05, 1.010810e-05, 1.003533e-05, 1.094816e-05, 1.174882e-05, 1.188400e-05, 1.206680e-05, 1.212674e-05, 1.249231e-05, 1.365744e-05, 1.632029e-05, 1.449828e-05, 1.629016e-05, 1.599958e-05, 1.968617e-05, 1.969068e-05, 1.972561e-05, 2.045977e-05, 1.836446e-05, 1.965865e-05, 2.016384e-05, 2.039124e-05, 2.408885e-05, 2.376285e-05, 2.395337e-05, 2.602860e-05, 2.862645e-05, 3.040946e-05, 3.216439e-05, 3.546673e-05, 3.527034e-05, 2.491410e-05, 2.462744e-05, 2.705848e-05, 2.722172e-05, 2.599691e-05, 2.776243e-05, 2.551111e-05, 2.705106e-05, 2.673012e-05, 2.795282e-05, 2.698715e-05, 3.002982e-05, 2.824995e-05, 2.534688e-05, 2.714108e-05, 2.850203e-05, 2.965501e-05, 3.131327e-05, 2.947131e-05, 3.253193e-05, 3.416903e-05, 3.058021e-05, 2.948300e-05, 3.239257e-05, 2.931725e-05, 3.484999e-05];

const stdday = [3.366163e-07, 7.089606e-07, 1.002567e-06, 1.467257e-06, 1.965123e-06, 2.079164e-06, 2.680844e-06, 3.431232e-06, 3.765651e-06, 3.850529e-06, 4.297061e-06, 5.300426e-06, 5.049050e-06, 6.632460e-06, 7.025610e-06, 8.815179e-06, 8.059801e-06, 9.133669e-06, 9.586492e-06, 1.026125e-05, 1.104400e-05, 1.280408e-05, 1.284779e-05, 1.338951e-05, 1.380421e-05, 1.605975e-05, 1.963515e-05, 2.235794e-05, 2.254865e-05, 2.664295e-05, 2.801855e-05, 2.702516e-05, 2.964451e-05, 3.364483e-05, 3.360693e-05, 3.651717e-05, 3.867344e-05, 3.882180e-05, 4.366121e-05, 4.877622e-05, 5.426914e-05, 6.014825e-05, 5.988709e-05, 6.388541e-05, 6.582952e-05, 7.080398e-05, 7.623408e-05, 7.978658e-05, 7.887853e-05, 8.872322e-05, 9.724181e-05, 9.071281e-05, 1.014742e-04, 7.900412e-05, 1.391666e-05, 1.406271e-05, 1.501573e-05, 1.501232e-05, 1.487840e-05, 1.503052e-05, 1.497957e-05, 1.505239e-05, 1.581055e-05, 1.584377e-05, 1.487225e-05, 1.499218e-05, 1.494824e-05, 1.502574e-05, 1.489518e-05, 1.499861e-05, 1.503177e-05, 1.483260e-05, 1.480940e-05, 1.540227e-05, 1.501651e-05, 1.490877e-05, 1.500109e-05, 1.503717e-05, 1.497237e-05, 1.497219e-05];

const ndrday = [2.444860e-07, 5.222806e-07, 7.442844e-07, 1.006310e-06, 1.480362e-06, 1.790019e-06, 2.601614e-06, 2.920642e-06, 3.497881e-06, 3.556354e-06, 3.835438e-06, 4.071524e-06, 4.704491e-06, 5.213196e-06, 5.633990e-06, 6.581441e-06, 7.230056e-06, 7.696886e-06, 7.894681e-06, 8.168179e-06, 8.471955e-06, 8.815676e-06, 9.268349e-06, 9.463110e-06, 1.067001e-05, 1.014601e-05, 1.094816e-05, 1.174882e-05, 1.183436e-05, 1.222415e-05, 1.216084e-05, 1.243616e-05, 1.365744e-05, 1.532009e-05, 1.456623e-05, 1.490676e-05, 1.553647e-05, 1.968617e-05, 1.652085e-05, 1.953990e-05, 2.045977e-05, 1.965865e-05, 1.997219e-05, 2.092483e-05, 2.054565e-05, 2.211872e-05, 2.269355e-05, 2.617345e-05, 2.561748e-05, 2.827255e-05, 2.825353e-05, 2.958450e-05, 2.775005e-05, 2.778078e-05, 2.998339e-05, 2.657015e-05, 2.929425e-05, 3.025374e-05, 2.665705e-05, 2.581609e-05, 2.665862e-05, 2.931089e-05, 2.609463e-05, 2.631122e-05, 2.964625e-05, 2.748686e-05, 2.532800e-05, 2.814977e-05, 2.666315e-05, 2.586391e-05, 2.697592e-05, 2.823705e-05, 2.758130e-05, 2.577017e-05, 2.676364e-05, 2.594511e-05, 2.652901e-05, 2.614594e-05, 2.674509e-05, 2.535949e-05];

const stdsmall = [3.744530e-07, 7.904913e-07, 1.300610e-06, 1.877682e-06, 2.359561e-06, 2.759019e-06, 3.542282e-06, 4.578075e-06, 4.956325e-06, 5.059711e-06, 5.441521e-06, 6.248918e-06, 6.373600e-06, 7.055431e-06, 8.704836e-06, 8.791013e-06, 9.453669e-06, 1.193246e-05, 1.115746e-05, 1.237196e-05, 1.288530e-05, 1.440082e-05, 1.459674e-05, 1.572306e-05, 1.746393e-05, 1.736392e-05, 2.213893e-05, 2.397940e-05, 2.327979e-05, 2.534251e-05, 2.756687e-05, 2.911142e-05, 3.103693e-05, 3.761810e-05, 3.774134e-05, 3.966310e-05, 4.710871e-05, 4.383507e-05, 4.435849e-05, 4.886874e-05, 5.722204e-05, 6.278839e-05, 6.479634e-05, 6.719520e-05, 7.248732e-05, 7.878058e-05, 7.797551e-05, 9.409875e-05, 9.399809e-05, 1.056246e-04, 1.011416e-04, 1.112904e-04, 1.253134e-04, 1.163310e-04, 1.397180e-04, 1.412168e-04, 1.634293e-04, 1.566923e-04, 1.864070e-04, 1.861402e-04, 1.874116e-04, 2.090663e-04, 2.353684e-04, 2.449997e-04, 2.666526e-04, 2.791808e-04, 3.244168e-04, 3.184042e-04, 3.630641e-04, 4.285631e-04, 4.673357e-04, 5.025180e-04, 6.348787e-04, 7.016320e-04, 7.809442e-04, 1.056494e-03, 1.324302e-03, 1.902986e-03, 3.815709e-03, 1.800000e+02];

const ndrsmall = [ 3.744530e-07, 7.849354e-07, 1.012175e-06, 1.487144e-06, 1.947993e-06, 2.078544e-06, 2.735251e-06, 3.504150e-06, 3.765651e-06, 3.912046e-06, 4.731274e-06, 4.713891e-06, 5.039958e-06, 6.629545e-06, 7.025610e-06, 6.907259e-06, 7.527196e-06, 7.769408e-06, 8.070465e-06, 8.542728e-06, 8.904458e-06, 8.888280e-06, 9.302651e-06, 1.037746e-05, 1.067001e-05, 1.073044e-05, 1.122700e-05, 1.180304e-05, 1.304042e-05, 1.315478e-05, 1.547164e-05, 1.542429e-05, 1.593398e-05, 1.637267e-05, 1.824440e-05, 1.781887e-05, 1.951870e-05, 1.968617e-05, 2.000454e-05, 2.073200e-05, 2.048995e-05, 2.173266e-05, 2.065080e-05, 2.239318e-05, 2.527703e-05, 2.456203e-05, 2.610315e-05, 2.854015e-05, 2.879501e-05, 3.323136e-05, 3.246935e-05, 3.681896e-05, 3.546673e-05, 3.436938e-05, 4.071821e-05, 4.331668e-05, 4.321151e-05, 4.483284e-05, 5.395275e-05, 5.635734e-05, 5.521645e-05, 5.666057e-05, 6.232493e-05, 6.371053e-05, 6.859644e-05, 8.110929e-05, 7.962349e-05, 9.539721e-05, 9.255844e-05, 1.077339e-04, 1.262085e-04, 1.356423e-04, 1.861125e-04, 1.724116e-04, 2.336238e-04, 2.803281e-04, 3.430796e-04, 5.399669e-04, 1.053252e-03, 1.800000e+02];

const ndrbf0 = [8.105798e-03, 8.740533e-03, 1.924753e-02, 1.798588e-02, 1.737519e-02, 1.155669e-02, 1.362648e-02, 2.135117e-02, 1.673397e-02, 1.838899e-02, 1.720053e-02, 1.928847e-02, 1.644710e-02, 1.948183e-02, 1.887864e-02, 2.344121e-02, 1.755570e-02, 2.144679e-02, 2.035396e-02, 2.106430e-02, 2.331826e-02, 1.972770e-02, 2.159706e-02, 2.073645e-02, 2.506679e-02, 2.107796e-02, 2.376322e-02, 2.323793e-02, 1.988951e-02, 2.255328e-02, 2.753932e-02, 2.394665e-02, 2.753932e-02, 2.277200e-02, 2.737541e-02, 2.770325e-02, 2.699293e-02, 2.697933e-02, 2.605034e-02, 2.509412e-02, 2.764498e-02, 2.599949e-02, 3.194192e-02, 2.793547e-02, 2.901805e-02, 2.862852e-02, 3.684766e-02, 2.967208e-02, 3.257013e-02, 3.328048e-02, 3.255598e-02, 3.347431e-02, 3.187583e-02, 3.785902e-02, 3.885350e-02, 3.840505e-02, 2.986538e-02, 3.404546e-02, 3.377204e-02, 3.443369e-02, 3.228332e-02, 3.694107e-02, 3.532896e-02, 3.355949e-02, 3.406007e-02, 3.951694e-02, 3.415706e-02, 3.756846e-02, 3.338976e-02, 3.560240e-02, 3.747452e-02, 3.542479e-02, 3.375871e-02, 3.665422e-02, 3.874157e-02, 3.381750e-02, 3.593030e-02, 3.929174e-02, 3.415763e-02, 4.777324e-02];

const stdbf0 = [2.804166e-03, 4.417325e-03, 5.712151e-03, 6.486021e-03, 1.017214e-02, 1.214567e-02, 1.455041e-02, 1.259487e-02, 8.814446e-03, 1.175275e-02, 1.316863e-02, 1.391992e-02, 1.748739e-02, 1.971404e-02, 1.449577e-02, 1.870317e-02, 1.974136e-02, 1.353743e-02, 1.449577e-02, 1.900370e-02, 1.924958e-02, 1.737812e-02, 1.939985e-02, 1.899004e-02, 1.730981e-02, 1.810211e-02, 1.610770e-02, 1.908566e-02, 1.811577e-02, 1.929057e-02, 1.909932e-02, 1.952279e-02, 1.823872e-02, 1.810211e-02, 1.942717e-02, 1.967306e-02, 1.879879e-02, 1.975502e-02, 1.941351e-02, 1.953645e-02, 1.375600e-02, 1.396091e-02, 1.386528e-02, 1.379699e-02, 1.375601e-02, 1.375602e-02, 1.375600e-02, 1.334619e-02, 1.338717e-02, 1.371502e-02, 1.161133e-02, 1.043655e-02, 1.310486e-02, 1.244461e-02, 9.958426e-03, 9.878971e-03, 1.003520e-02, 1.015814e-02, 1.036748e-02, 1.019912e-02, 1.027091e-02, 1.095271e-02, 1.299902e-02, 1.092815e-02, 9.788404e-03, 1.018546e-02, 1.016334e-02, 1.013608e-02, 1.149155e-02, 9.903852e-03, 9.871321e-03, 9.885201e-03, 1.007638e-02, 1.031387e-02, 1.004899e-02, 1.151366e-02, 1.141915e-02, 1.038773e-02, 1.007002e-02, 8.220133e-02];

const ndrok = [ 3.314157e-07, 7.414798e-07, 1.340293e-06, 1.500336e-06, 2.347513e-06, 2.687789e-06, 3.523046e-06, 4.355022e-06, 4.544263e-06, 4.665053e-06, 5.187013e-06, 5.642665e-06, 5.648806e-06, 6.923999e-06, 8.384703e-06, 7.804018e-06, 8.071463e-06, 8.228471e-06, 8.833536e-06, 8.943822e-06, 9.751226e-06, 1.010155e-05, 1.016610e-05, 1.143110e-05, 1.219638e-05, 1.082219e-05, 1.105898e-05, 1.291928e-05, 1.305886e-05, 1.308372e-05, 1.360181e-05, 1.410747e-05, 1.446881e-05, 1.439773e-05, 1.451214e-05, 1.481149e-05, 1.484694e-05, 1.497899e-05, 1.490520e-05, 1.469102e-05, 1.448845e-05,1.443887e-05, 1.450755e-05, 1.457813e-05, 1.446025e-05, 1.445342e-05, 1.509842e-05, 1.494390e-05, 1.502035e-05, 1.506726e-05, 1.490563e-05, 1.482591e-05, 1.459507e-05, 1.438928e-05, 1.359180e-05, 1.399279e-05, 1.398051e-05, 1.395092e-05, 1.406243e-05, 1.489813e-05, 1.404952e-05, 1.453042e-05, 1.481668e-05, 1.483929e-05, 1.429275e-05, 1.472702e-05, 1.476810e-05, 1.479279e-05, 1.441032e-05, 1.427637e-05, 1.428810e-05, 1.470500e-05, 1.427477e-05, 1.439435e-05, 1.441298e-05, 3.430035e-05, 2.895978e-05, 5.462125e-05, 1.092997e-04, 1.775642e+02];

const ndrokc = [2.622678e-07, 6.918605e-07, 1.167489e-06, 1.482148e-06, 2.347513e-06, 2.397669e-06, 2.671694e-06, 3.688556e-06, 4.544263e-06, 4.654336e-06, 5.187013e-06, 5.029108e-06, 5.646571e-06, 6.887813e-06, 7.116506e-06, 7.804018e-06, 7.800468e-06, 8.180182e-06, 8.347558e-06, 8.891348e-06, 9.751226e-06, 9.178066e-06, 1.014320e-05, 9.539234e-06, 1.058258e-05, 9.965329e-06, 1.094116e-05, 1.147091e-05, 1.297188e-05, 1.065140e-05, 1.360181e-05, 1.355516e-05, 1.383549e-05, 1.403085e-05, 1.446595e-05, 1.418174e-05, 1.468797e-05, 1.497899e-05, 1.443395e-05, 1.443583e-05, 1.133213e-05,1.443887e-05, 1.450504e-05, 1.457813e-05, 1.446025e-05, 1.436301e-05, 1.427535e-05, 1.474006e-05, 1.483820e-05, 1.506726e-05, 1.464220e-05, 1.482591e-05, 1.431196e-05, 1.438928e-05, 1.297184e-05, 1.396755e-05, 1.398051e-05, 1.395092e-05, 1.406243e-05, 1.371549e-05, 1.385529e-05, 1.383817e-05, 1.481668e-05, 1.412125e-05, 1.412862e-05, 1.443039e-05, 1.476810e-05, 1.423124e-05, 1.425124e-05, 1.415157e-05, 1.398697e-05, 1.435195e-05, 1.427477e-05, 1.425102e-05, 1.441298e-05, 1.436780e-05, 1.437925e-05, 1.416659e-05, 1.424719e-05, 1.412792e-05];

const stdokc = [3.346917e-07, 9.578924e-07, 1.210625e-06, 1.874781e-06, 2.505848e-06, 2.746285e-06, 3.565667e-06, 4.198671e-06, 4.207518e-06, 5.892483e-06, 5.166599e-06, 6.178459e-06, 5.836336e-06, 6.952337e-06, 7.167504e-06, 7.906112e-06, 8.445097e-06, 8.918862e-06, 8.901366e-06, 1.153221e-05, 1.155563e-05, 1.189788e-05, 1.224067e-05, 1.214075e-05, 1.271677e-05, 1.329759e-05, 1.334520e-05, 1.373724e-05, 1.884878e-05, 1.528106e-05, 1.915658e-05, 1.751454e-05, 1.600944e-05, 1.946288e-05, 2.046274e-05, 1.962320e-05, 1.987813e-05, 1.965406e-05, 2.050541e-05, 2.048111e-05, 2.001425e-05, 2.069134e-05, 2.386750e-05, 2.396392e-05, 2.390969e-05, 2.383971e-05, 2.372293e-05, 2.390427e-05, 2.360719e-05, 2.360852e-05, 2.329091e-05, 2.344288e-05, 2.590365e-05, 2.297921e-05, 2.605044e-05, 2.193215e-05, 2.212392e-05, 2.126918e-05, 2.218603e-05, 2.092305e-05, 2.023396e-05, 1.959073e-05, 1.930593e-05, 1.892217e-05, 1.819799e-05, 1.647049e-05, 1.598383e-05, 1.608193e-05, 1.476274e-05, 1.409724e-05, 1.444056e-05, 1.418945e-05, 1.453471e-05, 1.434607e-05, 1.420645e-05, 1.435723e-05, 1.401596e-05, 1.433243e-05, 1.446528e-05, 2.220616e-05];

const ndrokcd = [2.150735e-07, 4.798970e-07, 7.344784e-07, 9.751507e-07, 1.519374e-06, 1.810196e-06, 2.568172e-06, 2.715876e-06, 2.987993e-06, 3.182016e-06, 3.851741e-06, 3.916360e-06, 4.071319e-06, 5.118302e-06, 5.341728e-06, 5.715369e-06, 5.890366e-06, 6.281735e-06, 6.347660e-06, 6.463920e-06, 7.089652e-06, 7.206331e-06, 7.567520e-06, 7.652937e-06, 7.711626e-06, 7.770182e-06, 7.844697e-06, 8.593751e-06, 8.639733e-06, 8.670118e-06, 8.095813e-06, 8.138481e-06, 7.999749e-06, 7.775155e-06, 7.787467e-06, 7.922329e-06, 7.644994e-06, 7.722453e-06, 7.281120e-06, 7.949267e-06, 6.743197e-06,7.152271e-06, 6.798513e-06, 7.058266e-06, 6.911557e-06, 7.291439e-06, 7.088237e-06, 7.782632e-06, 7.044103e-06, 7.141886e-06, 7.131240e-06, 7.746183e-06, 7.598796e-06, 7.837137e-06, 7.929945e-06, 7.971435e-06, 8.124628e-06, 8.047090e-06, 7.965800e-06, 8.116940e-06, 7.966997e-06, 8.209299e-06, 8.808540e-06, 8.672490e-06, 8.801275e-06, 8.182235e-06, 8.013931e-06, 9.567472e-06, 8.104578e-06, 8.094826e-06, 8.013952e-06, 8.084668e-06, 8.067947e-06, 8.078949e-06, 8.113447e-06, 7.983502e-06, 7.998848e-06, 8.011594e-06, 7.998239e-06, 7.855958e-06];

const stdokcd = [2.150735e-07, 5.826456e-07, 1.010033e-06, 1.482497e-06, 1.604715e-06, 2.048760e-06, 2.647408e-06, 2.946316e-06, 3.112935e-06, 3.368460e-06, 3.827638e-06, 3.963443e-06, 4.189365e-06, 5.222807e-06, 5.619278e-06, 5.760912e-06, 5.917723e-06, 6.285933e-06, 6.949478e-06, 7.324895e-06, 7.423741e-06, 7.507527e-06, 7.588039e-06, 8.357825e-06, 8.422586e-06, 9.506772e-06, 9.503051e-06, 9.721047e-06, 1.311408e-05, 1.322163e-05, 1.369535e-05, 1.371002e-05, 1.079181e-05, 1.442962e-05, 1.468442e-05, 1.473637e-05, 1.487859e-05, 1.641870e-05, 1.647055e-05, 1.946899e-05, 1.933737e-05, 1.967526e-05, 2.040754e-05, 2.055405e-05, 2.035740e-05, 2.043611e-05, 2.147282e-05, 2.332474e-05, 2.120911e-05, 2.315832e-05, 2.286725e-05, 2.308380e-05, 2.257517e-05, 2.155470e-05, 1.951293e-05, 1.931687e-05, 1.968647e-05, 1.875354e-05, 1.878693e-05, 1.865996e-05, 1.787704e-05, 1.790868e-05, 1.817779e-05, 1.773547e-05, 1.634743e-05, 1.625142e-05, 1.401446e-05, 1.347582e-05, 1.341843e-05, 1.280223e-05, 1.216342e-05, 1.104637e-05, 1.033817e-05, 1.013448e-05, 9.923519e-06, 9.954146e-06, 9.690476e-06, 9.758583e-06, 9.808945e-06, 1.551398e-05];										    

const ndrdokcd = [2.134136e-07, 4.397891e-07, 4.769353e-07, 8.779406e-07, 9.537191e-07, 1.200415e-06, 1.696365e-06, 1.751541e-06, 1.903370e-06, 2.414833e-06, 2.414804e-06, 2.364633e-06, 2.414069e-06, 3.317194e-06, 3.424209e-06, 3.636982e-06, 3.753085e-06, 3.762884e-06, 3.754328e-06, 3.760613e-06, 4.452546e-06, 4.409318e-06, 4.425895e-06, 4.834238e-06, 4.848319e-06, 4.516977e-06, 5.131008e-06, 4.829650e-06, 4.825369e-06, 5.802236e-06, 5.880484e-06, 5.877874e-06, 6.288186e-06, 5.902796e-06, 5.920645e-06, 5.892178e-06, 5.908804e-06, 5.961328e-06, 6.114804e-06, 5.954242e-06, 5.968258e-06, 5.926593e-06, 6.004008e-06, 5.981887e-06, 6.324492e-06, 6.245957e-06, 6.165852e-06, 6.239144e-06, 6.216811e-06, 6.331209e-06, 6.364432e-06, 6.520594e-06, 6.460544e-06, 6.388309e-06, 5.939435e-06, 6.242718e-06, 6.184727e-06, 5.946396e-06, 5.928545e-06, 5.931382e-06, 5.913528e-06, 6.027904e-06, 6.161718e-06, 6.309155e-06, 6.208742e-06, 6.301214e-06, 6.437354e-06, 6.154315e-06, 6.114268e-06, 6.159352e-06, 6.319638e-06, 6.198569e-06, 6.085949e-06, 6.128577e-06, 6.145197e-06, 6.246214e-06, 6.193231e-06, 6.361269e-06, 6.283507e-06, 6.531273e-06];

const plotw = 900;
const ploth = 500;

const trace0 = {y:stdstd,   x:bax, name:'std/std' };
const trace1 = {y:stdday,   x:bax, name:'std/day' };
const trace2 = {y:stdsmall, x:bax, name:'std/small'};
const trace3 = {y:ndrstd,   x:bax, name:'ndr/std' };
const trace4 = {y:ndrday,   x:bax, name:'ndr/day' };
const trace5 = {y:ndrsmall, x:bax, name:'ndr/small'};

const trace6 = {y:ndrok,    x:bax, name:'ndr/overkill(bf)'};
const trace7 = {y:ndrokc,   x:bax, name:'ndr/overkill'};
const trace8 = {y:stdokc,   x:bax, name:'std/overkill'};
const trace9 = {y:ndrokcd,  x:bax, name:'ndr/overkill(dp)'};
const trace10= {y:stdokcd,  x:bax, name:'std/overkill(dp)'};
const trace11= {y:ndrdokcd, x:bax, name:'ndr(dp)/overkill(dp)'};

const tracex1 = {y:stdbf0, x:bax, name:'std/bf0'};
const tracex2 = {y:ndrbf0, x:bax, name:'ndr/bf0'};

const traces  = [trace0,trace1,trace2,trace3,trace4,trace5];
const traces2 = [tracex1, tracex2];
const traces3 = [trace2, trace5];

const traces4 = [trace6, trace4];
const traces5 = [trace7, trace8, trace9, trace10, trace4,trace1];
const traces6 = [trace9,trace11];

const layout = {
  title:  'roundtrip angle error by rot angle',
  yaxis: { range: [0.0, 0.00013], hoverformat: 'g', exponentformat: 'power' },
  xaxis: { nticks: 20 },
  height: ploth,
  width:  plotw
};

const layout2 = {
  title:  'roundtrip angle error by rot angle',
  xaxis: { nticks: 20 },
  height: ploth,
  width:  plotw
};

const layout3 = {
  title:  'roundtrip angle error by rot angle',
  yaxis: { range: [0.0, 0.003], hoverformat: 'g', exponentformat: 'power' },
  xaxis: { nticks: 20 },
  height: ploth,
  width:  plotw
};

const layout4 = {
  title:  'roundtrip angle error by rot angle',
  yaxis: { hoverformat: 'g', exponentformat: 'power' },
  xaxis: { nticks: 20 },
  height: ploth,
  width:  plotw
};


Plotly.newPlot('fig1',  traces,  layout,  {displaylogo: false, autosizable: true});
Plotly.newPlot('fig2',  traces2, layout2, {displaylogo: false, autosizable: true});
Plotly.newPlot('fig3',  traces3, layout3, {displaylogo: false, autosizable: true});
Plotly.newPlot('fig4',  traces4, layout,  {displaylogo: false, autosizable: true});
Plotly.newPlot('fig5',  traces5, layout4, {displaylogo: false, autosizable: true});
Plotly.newPlot('fig7',  traces6, layout4, {displaylogo: false, autosizable: true});


const rate = {
x:[2.04082e-8, 0.000306718, 0.000613415, 0.00122681, 0.0024536, 0.00490718, 0.00981434, 0.0196287, 0.0409084, 0.0607779, 0.0802576, 0.101388, 0.121109, 0.142481, 0.163463, 0.183035, 0.204257, 0.22407, 0.243493, 0.264567, 0.284231, 0.305546, 0.326471, 0.345986, 0.367151, 0.386907, 0.408314, 0.429331, 0.448938, 0.470196, 0.490044, 0.509502, 0.530611, 0.55031, 0.57166, 0.59262, 0.61217, 0.633371, 0.653162, 0.672563, 0.693616, 0.713258, 0.734551, 0.754434, 0.773927, 0.795071, 0.814805, 0.83619, 0.857186, 0.876771, 0.898007, 0.917833, 0.93727, 0.937599, 0.937929, 0.938588, 0.939906, 0.942542, 0.947813, 0.948143, 0.948472, 0.949131, 0.950449, 0.953085, 0.958357, 0.958665, 0.958972, 0.959587, 0.960817, 0.963276, 0.968196, 0.968503, 0.968811, 0.969426, 0.970655, 0.973115, 0.978034, 0.978378, 0.978721, 0.979407, 0.98078, 0.983526, 0.983869, 0.984212, 0.984899, 0.986271, 0.989017, 0.98936, 0.989704, 0.99039, 0.991763, 0.992106, 0.992449, 0.993136, 0.994509, 0.994852, 0.995195, 0.995881, 0.996225, 0.996568, 0.997254, 0.997597, 0.997941, 0.998284, 0.998627, 0.99897, 0.999314, 0.999657, 1.],
y:[1., 0.999609, 0.999219, 0.998438, 0.996876, 0.993752, 0.987504, 0.97501, 0.947928, 0.922663, 0.897923, 0.87113, 0.846177, 0.819203, 0.792803, 0.768261, 0.741751, 0.717111, 0.693067, 0.667115, 0.64304, 0.617108, 0.591832, 0.568431, 0.543254, 0.519957, 0.494949, 0.47065, 0.448222, 0.424184, 0.402016, 0.380558, 0.357603, 0.336501, 0.313998, 0.2923, 0.272432, 0.251313, 0.232021, 0.213526, 0.193951, 0.176176, 0.157474, 0.140574, 0.124567, 0.107874, 0.0929669, 0.077602, 0.0633813, 0.0509575, 0.0384972, 0.0279223, 0.0186819, 0.0185358, 0.0183901, 0.0180999, 0.0175239, 0.0163903, 0.0141985, 0.014065, 0.0139318, 0.0136668, 0.0131418, 0.0121121, 0.0101371, 0.0100255, 0.00991433, 0.00969318, 0.00925583, 0.00840127, 0.0067761, 0.00667839, 0.00658115, 0.00638808, 0.00600767, 0.00527035, 0.00389507, 0.00380433, 0.00371431, 0.00353641, 0.00318939, 0.00253202, 0.00245344, 0.00237568, 0.00222267, 0.00192697, 0.0013794, 0.00131532, 0.00125225, 0.00112925, 0.000896317, 0.000840929, 0.000786727, 0.000681993, 0.000488096, 0.000443082, 0.000399542, 0.000317093, 0.000278309, 0.000241246, 0.00017264, 0.000141311, 0.000112145, 0.0000853162, 0.0000610507, 0.0000396559, 0.0000215874, 7.63303e-6, 3.49977e-12],
mode: 'lines'
};

Plotly.newPlot('fig6', [rate], {height: 400, width: 600}, {displaylogo: false, autosizable: true});


</script>
