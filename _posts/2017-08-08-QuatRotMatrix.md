---
layout:       post
title:        On quaternion/rotation matrix conversions and errors
categories:   [quaternions]
tags:         [rotations, conversion]
plotly:       true
disqus:       false
description:  Skims the math of quaternions to rotations matrices and back with empirical errors including small angle special case.
---

Toy code can be found: [HERE](http://github.com/Marc-B-Reynolds/Stand-alone-junk/blob/master/src/Posts/q2mat.c)
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
\mathbf{x} = \left(1,0,0\right) \\
\mathbf{y} = \left(0,1,0\right) \\
\mathbf{z} = \left(0,0,1\right) \\
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

Another aside: In cases where input quaternions are moving away from unity due to compounding errors it is possible to fold the corrective step into the transform instead of renormalizing as a pre-step.  As an example:

{% highlight c %}
void quat_to_mat33_nu(mat33_t* m, quat_t* q)
{
  float x  = q->x, y  = q->y, z  = q->z, w  = q->w;
  float xx = x*x,  yy = y*y,  zz = z*z,  ww = w*w;
  float s  = 2.f*recipf(xx+yy+zz+ww);
  float sx = s*x,  sy = s*y,  sz = s*z;
  float xy = sy*x, xz = sz*x, yz = sy*z;
  float wx = sx*w, wy = sy*w, wz = sz*w;

  m->m00 = 1.f-s*(yy+zz); m->m11 = 1.f-s*(xx+zz); m->m22 = 1.f-s*(xx+yy);

  m->m10 = xy+wz; m->m01 = xy-wz;
  m->m20 = xz-wy; m->m02 = xz+wy;
  m->m21 = yz+wx; m->m12 = yz-wx;
}
{% endhighlight %}

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
m_{00} < 0 \implies y^2+z^2 > \frac{1}{2} \\
m_{11} < 0 \implies x^2+z^2 > \frac{1}{2} \\
m_{22} < 0 \implies x^2+y^2 > \frac{1}{2}
$$

\\
and the permutations of sum/difference of pairs of diagonal elements imply:

$$
m_{00}+m_{11} < 0 \implies w^2 < z^2 \\
m_{00}-m_{11} < 0 \implies x^2 < y^2 \\
m_{00}+m_{22} < 0 \implies w^2 < y^2 \\
m_{00}-m_{22} < 0 \implies x^2 < z^2 \\ 
m_{11}+m_{22} < 0 \implies w^2 < x^2 \\
m_{11}-m_{22} < 0 \implies y^2 < z^2
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

<div id="fig1" style="width:100%"></div>

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

References and Footnotes
------

[^day]: *"Converting a Rotation Matrix to a Quaternion"*, Mike Day, 2015 ([download page](http://www.insomniacgames.com/converting-a-rotation-matrix-to-a-quaternion/))

<script>

var bax = new Array(80);

{ for(var i=0; i<80; i++) bax[i] = i/79.0; }

var stdstd = [3.366163e-07, 7.089606e-07, 1.002567e-06, 1.467112e-06, 1.965123e-06, 2.060151e-06, 2.680844e-06, 3.431232e-06, 3.765651e-06, 3.850529e-06, 4.297061e-06, 5.300426e-06, 5.037166e-06, 6.900605e-06, 7.025214e-06, 8.815179e-06, 8.059801e-06, 9.133669e-06, 1.013469e-05, 1.025556e-05, 1.104310e-05, 1.280408e-05, 1.284779e-05, 1.341579e-05, 1.380421e-05, 1.606909e-05, 1.868217e-05, 2.235794e-05, 2.254608e-05, 2.664295e-05, 2.799844e-05, 2.702300e-05, 2.964451e-05, 3.336401e-05, 3.319995e-05, 3.651717e-05, 3.869078e-05, 3.879974e-05, 4.364011e-05, 4.584912e-05, 5.431133e-05, 6.020936e-05, 6.017704e-05, 6.368387e-05, 6.587051e-05, 7.539200e-05, 8.067569e-05, 8.180513e-05, 8.839987e-05, 9.791440e-05, 1.004588e-04, 1.025081e-04, 1.215223e-04, 1.095087e-04, 1.541106e-05, 1.501028e-05, 1.511791e-05, 1.578418e-05, 1.541071e-05, 1.514047e-05, 1.575093e-05, 1.585129e-05, 1.577171e-05, 1.574984e-05, 1.596803e-05, 1.575143e-05, 1.604615e-05, 1.594718e-05, 1.552819e-05, 1.609223e-05, 1.552909e-05, 1.577294e-05, 1.581677e-05, 1.519899e-05, 1.591611e-05, 1.574458e-05, 1.565862e-05, 1.579813e-05, 1.572727e-05, 1.575847e-05];

var ndrstd = [2.444860e-07, 5.222806e-07, 7.442844e-07, 1.006310e-06, 1.459467e-06, 1.721504e-06, 2.588514e-06, 2.920642e-06, 3.033302e-06, 3.556354e-06, 3.844039e-06, 4.009854e-06, 4.704491e-06, 5.158537e-06, 5.627595e-06, 6.705708e-06, 7.230056e-06, 7.696886e-06, 7.819148e-06, 7.650115e-06, 8.471955e-06, 8.790127e-06, 9.268349e-06, 1.038301e-05, 1.010810e-05, 1.003533e-05, 1.094816e-05, 1.174882e-05, 1.188400e-05, 1.206680e-05, 1.212674e-05, 1.249231e-05, 1.365744e-05, 1.632029e-05, 1.449828e-05, 1.629016e-05, 1.599958e-05, 1.968617e-05, 1.969068e-05, 1.972561e-05, 2.045977e-05, 1.836446e-05, 1.965865e-05, 2.016384e-05, 2.039124e-05, 2.408885e-05, 2.376285e-05, 2.395337e-05, 2.602860e-05, 2.862645e-05, 3.040946e-05, 3.216439e-05, 3.546673e-05, 3.527034e-05, 2.491410e-05, 2.462744e-05, 2.705848e-05, 2.722172e-05, 2.599691e-05, 2.776243e-05, 2.551111e-05, 2.705106e-05, 2.673012e-05, 2.795282e-05, 2.698715e-05, 3.002982e-05, 2.824995e-05, 2.534688e-05, 2.714108e-05, 2.850203e-05, 2.965501e-05, 3.131327e-05, 2.947131e-05, 3.253193e-05, 3.416903e-05, 3.058021e-05, 2.948300e-05, 3.239257e-05, 2.931725e-05, 3.484999e-05];

var stdday = [3.366163e-07, 7.089606e-07, 1.002567e-06, 1.467257e-06, 1.965123e-06, 2.079164e-06, 2.680844e-06, 3.431232e-06, 3.765651e-06, 3.850529e-06, 4.297061e-06, 5.300426e-06, 5.049050e-06, 6.632460e-06, 7.025610e-06, 8.815179e-06, 8.059801e-06, 9.133669e-06, 9.586492e-06, 1.026125e-05, 1.104400e-05, 1.280408e-05, 1.284779e-05, 1.338951e-05, 1.380421e-05, 1.605975e-05, 1.963515e-05, 2.235794e-05, 2.254865e-05, 2.664295e-05, 2.801855e-05, 2.702516e-05, 2.964451e-05, 3.364483e-05, 3.360693e-05, 3.651717e-05, 3.867344e-05, 3.882180e-05, 4.366121e-05, 4.877622e-05, 5.426914e-05, 6.014825e-05, 5.988709e-05, 6.388541e-05, 6.582952e-05, 7.080398e-05, 7.623408e-05, 7.978658e-05, 7.887853e-05, 8.872322e-05, 9.724181e-05, 9.071281e-05, 1.014742e-04, 7.900412e-05, 1.391666e-05, 1.406271e-05, 1.501573e-05, 1.501232e-05, 1.487840e-05, 1.503052e-05, 1.497957e-05, 1.505239e-05, 1.581055e-05, 1.584377e-05, 1.487225e-05, 1.499218e-05, 1.494824e-05, 1.502574e-05, 1.489518e-05, 1.499861e-05, 1.503177e-05, 1.483260e-05, 1.480940e-05, 1.540227e-05, 1.501651e-05, 1.490877e-05, 1.500109e-05, 1.503717e-05, 1.497237e-05, 1.497219e-05];

var ndrday = [2.444860e-07, 5.222806e-07, 7.442844e-07, 1.006310e-06, 1.480362e-06, 1.790019e-06, 2.601614e-06, 2.920642e-06, 3.497881e-06, 3.556354e-06, 3.835438e-06, 4.071524e-06, 4.704491e-06, 5.213196e-06, 5.633990e-06, 6.581441e-06, 7.230056e-06, 7.696886e-06, 7.894681e-06, 8.168179e-06, 8.471955e-06, 8.815676e-06, 9.268349e-06, 9.463110e-06, 1.067001e-05, 1.014601e-05, 1.094816e-05, 1.174882e-05, 1.183436e-05, 1.222415e-05, 1.216084e-05, 1.243616e-05, 1.365744e-05, 1.532009e-05, 1.456623e-05, 1.490676e-05, 1.553647e-05, 1.968617e-05, 1.652085e-05, 1.953990e-05, 2.045977e-05, 1.965865e-05, 1.997219e-05, 2.092483e-05, 2.054565e-05, 2.211872e-05, 2.269355e-05, 2.617345e-05, 2.561748e-05, 2.827255e-05, 2.825353e-05, 2.958450e-05, 2.775005e-05, 2.778078e-05, 2.998339e-05, 2.657015e-05, 2.929425e-05, 3.025374e-05, 2.665705e-05, 2.581609e-05, 2.665862e-05, 2.931089e-05, 2.609463e-05, 2.631122e-05, 2.964625e-05, 2.748686e-05, 2.532800e-05, 2.814977e-05, 2.666315e-05, 2.586391e-05, 2.697592e-05, 2.823705e-05, 2.758130e-05, 2.577017e-05, 2.676364e-05, 2.594511e-05, 2.652901e-05, 2.614594e-05, 2.674509e-05, 2.535949e-05];

var stdsmall = [3.744530e-07, 7.904913e-07, 1.300610e-06, 1.877682e-06, 2.359561e-06, 2.759019e-06, 3.542282e-06, 4.578075e-06, 4.956325e-06, 5.059711e-06, 5.441521e-06, 6.248918e-06, 6.373600e-06, 7.055431e-06, 8.704836e-06, 8.791013e-06, 9.453669e-06, 1.193246e-05, 1.115746e-05, 1.237196e-05, 1.288530e-05, 1.440082e-05, 1.459674e-05, 1.572306e-05, 1.746393e-05, 1.736392e-05, 2.213893e-05, 2.397940e-05, 2.327979e-05, 2.534251e-05, 2.756687e-05, 2.911142e-05, 3.103693e-05, 3.761810e-05, 3.774134e-05, 3.966310e-05, 4.710871e-05, 4.383507e-05, 4.435849e-05, 4.886874e-05, 5.722204e-05, 6.278839e-05, 6.479634e-05, 6.719520e-05, 7.248732e-05, 7.878058e-05, 7.797551e-05, 9.409875e-05, 9.399809e-05, 1.056246e-04, 1.011416e-04, 1.112904e-04, 1.253134e-04, 1.163310e-04, 1.397180e-04, 1.412168e-04, 1.634293e-04, 1.566923e-04, 1.864070e-04, 1.861402e-04, 1.874116e-04, 2.090663e-04, 2.353684e-04, 2.449997e-04, 2.666526e-04, 2.791808e-04, 3.244168e-04, 3.184042e-04, 3.630641e-04, 4.285631e-04, 4.673357e-04, 5.025180e-04, 6.348787e-04, 7.016320e-04, 7.809442e-04, 1.056494e-03, 1.324302e-03, 1.902986e-03, 3.815709e-03, 1.800000e+02];

var ndrsmall = [ 3.744530e-07, 7.849354e-07, 1.012175e-06, 1.487144e-06, 1.947993e-06, 2.078544e-06, 2.735251e-06, 3.504150e-06, 3.765651e-06, 3.912046e-06, 4.731274e-06, 4.713891e-06, 5.039958e-06, 6.629545e-06, 7.025610e-06, 6.907259e-06, 7.527196e-06, 7.769408e-06, 8.070465e-06, 8.542728e-06, 8.904458e-06, 8.888280e-06, 9.302651e-06, 1.037746e-05, 1.067001e-05, 1.073044e-05, 1.122700e-05, 1.180304e-05, 1.304042e-05, 1.315478e-05, 1.547164e-05, 1.542429e-05, 1.593398e-05, 1.637267e-05, 1.824440e-05, 1.781887e-05, 1.951870e-05, 1.968617e-05, 2.000454e-05, 2.073200e-05, 2.048995e-05, 2.173266e-05, 2.065080e-05, 2.239318e-05, 2.527703e-05, 2.456203e-05, 2.610315e-05, 2.854015e-05, 2.879501e-05, 3.323136e-05, 3.246935e-05, 3.681896e-05, 3.546673e-05, 3.436938e-05, 4.071821e-05, 4.331668e-05, 4.321151e-05, 4.483284e-05, 5.395275e-05, 5.635734e-05, 5.521645e-05, 5.666057e-05, 6.232493e-05, 6.371053e-05, 6.859644e-05, 8.110929e-05, 7.962349e-05, 9.539721e-05, 9.255844e-05, 1.077339e-04, 1.262085e-04, 1.356423e-04, 1.861125e-04, 1.724116e-04, 2.336238e-04, 2.803281e-04, 3.430796e-04, 5.399669e-04, 1.053252e-03, 1.800000e+02];

var ndrbf0 = [8.105798e-03, 8.740533e-03, 1.924753e-02, 1.798588e-02, 1.737519e-02, 1.155669e-02, 1.362648e-02, 2.135117e-02, 1.673397e-02, 1.838899e-02, 1.720053e-02, 1.928847e-02, 1.644710e-02, 1.948183e-02, 1.887864e-02, 2.344121e-02, 1.755570e-02, 2.144679e-02, 2.035396e-02, 2.106430e-02, 2.331826e-02, 1.972770e-02, 2.159706e-02, 2.073645e-02, 2.506679e-02, 2.107796e-02, 2.376322e-02, 2.323793e-02, 1.988951e-02, 2.255328e-02, 2.753932e-02, 2.394665e-02, 2.753932e-02, 2.277200e-02, 2.737541e-02, 2.770325e-02, 2.699293e-02, 2.697933e-02, 2.605034e-02, 2.509412e-02, 2.764498e-02, 2.599949e-02, 3.194192e-02, 2.793547e-02, 2.901805e-02, 2.862852e-02, 3.684766e-02, 2.967208e-02, 3.257013e-02, 3.328048e-02, 3.255598e-02, 3.347431e-02, 3.187583e-02, 3.785902e-02, 3.885350e-02, 3.840505e-02, 2.986538e-02, 3.404546e-02, 3.377204e-02, 3.443369e-02, 3.228332e-02, 3.694107e-02, 3.532896e-02, 3.355949e-02, 3.406007e-02, 3.951694e-02, 3.415706e-02, 3.756846e-02, 3.338976e-02, 3.560240e-02, 3.747452e-02, 3.542479e-02, 3.375871e-02, 3.665422e-02, 3.874157e-02, 3.381750e-02, 3.593030e-02, 3.929174e-02, 3.415763e-02, 4.777324e-02];

var stdbf0 = [2.804166e-03, 4.417325e-03, 5.712151e-03, 6.486021e-03, 1.017214e-02, 1.214567e-02, 1.455041e-02, 1.259487e-02, 8.814446e-03, 1.175275e-02, 1.316863e-02, 1.391992e-02, 1.748739e-02, 1.971404e-02, 1.449577e-02, 1.870317e-02, 1.974136e-02, 1.353743e-02, 1.449577e-02, 1.900370e-02, 1.924958e-02, 1.737812e-02, 1.939985e-02, 1.899004e-02, 1.730981e-02, 1.810211e-02, 1.610770e-02, 1.908566e-02, 1.811577e-02, 1.929057e-02, 1.909932e-02, 1.952279e-02, 1.823872e-02, 1.810211e-02, 1.942717e-02, 1.967306e-02, 1.879879e-02, 1.975502e-02, 1.941351e-02, 1.953645e-02, 1.375600e-02, 1.396091e-02, 1.386528e-02, 1.379699e-02, 1.375601e-02, 1.375602e-02, 1.375600e-02, 1.334619e-02, 1.338717e-02, 1.371502e-02, 1.161133e-02, 1.043655e-02, 1.310486e-02, 1.244461e-02, 9.958426e-03, 9.878971e-03, 1.003520e-02, 1.015814e-02, 1.036748e-02, 1.019912e-02, 1.027091e-02, 1.095271e-02, 1.299902e-02, 1.092815e-02, 9.788404e-03, 1.018546e-02, 1.016334e-02, 1.013608e-02, 1.149155e-02, 9.903852e-03, 9.871321e-03, 9.885201e-03, 1.007638e-02, 1.031387e-02, 1.004899e-02, 1.151366e-02, 1.141915e-02, 1.038773e-02, 1.007002e-02, 8.220133e-02];



var trace0 = {y:stdstd,   x:bax, name:'std/std' };
var trace1 = {y:stdday,   x:bax, name:'std/day' };
var trace2 = {y:stdsmall, x:bax, name:'std/small'};
var trace3 = {y:ndrstd,   x:bax, name:'ndr/std' };
var trace4 = {y:ndrday,   x:bax, name:'ndr/day' };
var trace5 = {y:ndrsmall, x:bax, name:'ndr/small'};
//var trace7 = {y:foo,  x:bax, name:'foo'};

var tracex1 = {y:stdbf0, x:bax, name:'std/bf0'};
var tracex2 = {y:ndrbf0, x:bax, name:'ndr/bf0'};

var traces  = [trace0,trace1,trace2,trace3,trace4,trace5];
var traces2 = [tracex1, tracex2];
var traces3 = [trace2, trace5];

var layout = {
  title:  'roundtrip angle error by rot angle',
  yaxis: { range: [0.0, 0.00013], hoverformat: 'g', exponentformat: 'power' },
  xaxis: { nticks: 20 },
  height: 500,
  width:  800
};

var layout2 = {
  xaxis: { nticks: 20 },
  height: 500,
  width:  800
};

var layout3 = {
  yaxis: { range: [0.0, 0.003], hoverformat: 'g', exponentformat: 'power' },
  xaxis: { nticks: 20 },
  height: 500,
  width:  800
};


Plotly.newPlot('fig1',  traces,  layout,  {displaylogo: false, autosizable: true});
Plotly.newPlot('fig2',  traces2, layout2, {displaylogo: false, autosizable: true});
Plotly.newPlot('fig3',  traces3, layout3, {displaylogo: false, autosizable: true});

</script>
