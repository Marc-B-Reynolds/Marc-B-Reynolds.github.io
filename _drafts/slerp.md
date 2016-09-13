---
layout:       post
title:        On certain results of quaternionic interpolation
tagline:      (or let's walk a 2D path)
categories:   [quaternions]
tags:         [interpolation]
description:  slerp is dead! long live slerp!
plotly:       true
---

This note is on basic quaternion interpolation with the following goals:

* Derive slerp, attempt to make it less black-box like and show it is a 2D problem
* Show that it is equivalent to the commonly seen function
* Form reference versions that only require one forward trig operation
* Present analysis tools which require only reals
* Give basic analysis for lerp
* Show some known and relatively unknown formulations

Non-goals are to convince anyone that any specific choice is the correct answer nor to present efficent implementations.

The following shadertoy is a viz for three interpolations formulated directly in 2D: slerp (blue), nlerp (magenta) and half-angle+nlerp (red).  Clicking down sets the initial point and the final point is the current mouse position.  On button release the initial position resets to zero and the final remains.

<iframe width="640" height="360" frameborder="0" src="https://www.shadertoy.com/embed/4tdGWB?gui=true&t=10&paused=false&muted=false" allowfullscreen></iframe>


<br>

------

From complex real powers <small>to the algebraic form of slerp</small>
------

\\
Since quaternions are an extension of complex numbers let's first consider the complex case.  Given a unit complex number $z$ both in explict and polar form:

$$
\begin{eqnarray}
  z & = & a + b\mathbf{i} \label{zr} \\ 
    & = & \cos\left(\theta\right)+\sin\left(\theta\right)~\mathbf{i} \label{zp}
\end{eqnarray}
$$

\\
The real (principle) power can be expressed in terms of the polar form as $\left( t\in \mathbb{R} \right) $:

$$ \begin{equation} \label{power}
z^t = \cos\left(t\theta\right)+\sin\left(t\theta\right)~\mathbf{i}
\end{equation} $$ 


\\
This linearly parameterizes the angle between one and $z$.  If $z$ respresents a rotation then this produces a continuous interpolation between the identity and the rotation represented by $z$ following the minimal magnitude path.

To generalize to an arbitrary starting point $A$ and ending point $B$ we can simply transform the problem:

$$ \begin{equation} \label{zslerp}
  s \left( t \right) = \left(BA^*\right)^tA
\end{equation} $$ 

Which is the algebraic form of slerp in two-dimensions.


<br>

------

SLERP <small>the dreaded small angle circular arc</small>
------

\\
Replacing $z$ with a unit quaternion:

$$
\begin{eqnarray}
  Q & = & a + \left(x,~y,~z\right) \label{qxform} \\
    & = & a + b~\mathbf{u} \label{qiform} \\
    & = & \cos\left(\theta\right)+\sin\left(\theta\right)~\mathbf{u} \label{qpform}
\end{eqnarray}
$$

changes hardly anything.  I'm assuming familarity with the contents of [half/double angle and Cayley: preliminaries]({{site.base}}/quaternions/2016/05/30/QuatHDAngleCayley.html#prelim), the bullet points of which are:

* effectively working with the complex number: $z = a + b~\mathbf{i} $
* unit magnitude: $ a^2 + b^2 =1 $
* positive scalar: $ a $ on $\left[0,1\right] $
* $ b $ on $\left[0,1\right] $, since $b=\sqrt{x^2+y^2+z^2}$. Any sign of $b$ is contained by the implied $\mathbf{u}$
* the implied angle $\theta$ is therefore on $\left[0,~\frac{\pi}{2}\right]$
* the domain is the arc of unit circle in the first quadrant, axes included.

Since the quaternion product doesn't commute we can express slerp as:

$$
\begin{eqnarray}
  s\left( t \right) & = & \left(BA^*\right)^tA     \label{slerp}  \\
                    & = & A\left(A^*B\right)^t     \label{slerp1} \\
		    & = & \left(AB^*\right)^{1-t}B \label{slerp2} \\
		    & = & B\left(B^*A\right)^{1-t} \label{slerp3} 
\end{eqnarray}
$$

\\
Where $\eqref{slerp2}$ and $\eqref{slerp3}$ are simply reversing the points and parameterization (holds in complex as well).  Using $\eqref{slerp}$ we have $Q=BA^*$ where $A \cdot B \geq 0$ and the action of slerp is $Q^t$, which maps a coordinate in the {1,$Q$}-plane to a coordinate in the {1,$Q$}-plane.

<br>

------

Algebraic to geometric form <small>the version you almost always see</small>
------

\\
A skippable sketch of a grind through the math (longer than needed by not having a table of operators and identities):

$$ \begin{align*}
A & = \cos(\alpha) + \sin(\alpha)\mathbf{v}      \\
B & = \cos(\beta)  + \sin(\beta)\mathbf{w}       \\
                                                 \\
Q &= BA^{*} \\
  & = \cos(\alpha)\cos(\beta)+ \sin(\alpha)\sin(\beta)\left( \mathbf{v} \cdot \mathbf{w} \right) +
  \sin(\alpha)\sin(\beta)(\mathbf{v} \times \mathbf{w}) -
  \sin(\alpha)\cos(\beta) \mathbf{v} +
  \cos(\alpha)\sin(\beta) \mathbf{w}
  \\
   &= \frac{1}{2}\left(BA^*+AB^* \right) ~+ \frac{1}{2}\left(BA^*-AB^* \right) \\
   &= A \cdot B ~+ \frac{1}{2}\left(BA^*-AB^* \right) \\
   & =  \cos(\theta) + \sin(\theta)\mathbf{u}    \\
  \\
Q^t &= \cos(t \theta) + \sin(t \theta)\mathbf{u}  \\
   \\  
  s\left( t \right) & = Q^t A \\
    &= \left( \cos(t \theta) + \sin(t \theta)\mathbf{u} \right)A \\
    &= \cos(t \theta)A + \sin(t \theta)\mathbf{u}A \\


\end{align*} $$

\\
Or we can directly place our 2D coordinate $\eqref{power}$ back into the original space:


$$ \begin{equation} \label{slerpxy}
s(t) = \cos\left(t\theta\right)X + \sin\left(t\theta\right)Y
\end{equation} $$ 

\\
Note that this holds for any number of dimensions. 

Using $\eqref{slerp}$ or $\eqref{slerp1}$ then the ${X}$ direction is $A$ and ${Y}$ is the direction orthogonal to $A$ that contains $B$ (normalized).  Since the inputs are normalized we have:

$$ 
Y = \left\Vert B - (A \cdot B)A\right\Vert
$$

\\
which brings us to:

$$ \begin{equation} \label{srecon}
s(t) = \cos\left(t\theta\right)A + \sin\left(t\theta\right)\frac{B - \left(A\cdot B\right)A}{\sqrt{1 -\left(A\cdot B\right)^{2}}}
\end{equation} $$ 

\\
If we swap $A$ and $B$ and replacing $t$ by $1-t$ to instead use $\eqref{slerp2}$ or $\eqref{slerp3}$.  Making some trig replacements:

$$
\begin{eqnarray*}
  \cos\left(\theta\right) & = & A\cdot B  \\
  \sin\left(\theta\right) & = & \sqrt{1 -\left(A\cdot B\right)^{2}}
\end{eqnarray*}
$$

\\
then we have:

$$
\begin{eqnarray}
s(t) & = & \cos(t\theta) A + \sin(t\theta) \frac{B - A \cos(\theta)}{\sin(\theta)}  \nonumber                \\
     & = & \left(\cos(t\theta) - \frac{\cos(\theta)\sin(t\theta)}{sin(\theta)}\right)A + \frac{\sin(t\theta)}{\sin(\theta)}B  \nonumber     \\
     & = & \frac{\cos(t\theta)\sin(\theta) - \cos(\theta)\sin(t\theta)}{\sin(\theta)}A + \frac{\sin(t\theta)}{\sin(\theta)}B \nonumber \\
     & = & \frac{\sin\left(\left(1-t\right)\theta\right)}{\sin\left(\theta\right)}A+\frac{\sin\left(t\theta\right)}{\sin\left(\theta\right)}B \label{gdavis}
\end{eqnarray}
$$

\\
Where $\eqref{gdavis}$ can be directly formulated from geometric reasoning and has somehow has be taken to be the *true* defination of slerp.  

Ken Shoemake's 1985 paper[^shoemake85] gives $\eqref{slerp1}$ and $\eqref{gdavis}$.  This equation has be covered very well, as an example J.M.P van Waveren[^id] 2005.

<br>

------

Reference formulation of slerp <small>one forward trig op is enough</small>
------

\\
Starting from $\eqref{srecon}$, replace $A\cdot B$ by $d$ and collecting terms gives:

$$
\begin{eqnarray*}
s(t) & = & \cos\left(t\theta\right)A + \sin\left(t\theta\right)\frac{B - dA}{\sqrt{1 -d^2}} \\
    & = & \left( \cos\left(t\theta\right)-d\frac{\sin\left(t\theta\right)}{\sqrt{1 -d^2}} \right) A + \frac{\sin\left(t\theta\right)}{\sqrt{1 -d^2}}B
\end{eqnarray*}
$$

\\
The two forward trig ops can be converted into one and given the range of both are positive[^nosign].  Showing both options as a pair of scale values where $s_t=\sin\left(t\theta\right)$ and $c_t=\cos\left(t\theta\right)$:

$$
\begin{eqnarray*}
s(t) & = & \left\{\sqrt{1-s_t^2}-d \frac{s_t}{\sqrt{1-d^2}},~\frac{s_t}{\sqrt{1-d^2}}\right\} \\
     & = & \left\{c_t-d \sqrt{\frac{1-c_t^2}{1-d^2}},~\sqrt{\frac{1-c_t^2}{1-d^2}}\right\}
\end{eqnarray*}
$$

\\
This leaves the computation of theta. I'll paritially dismiss (directly) using $\text{acos}$ and $\text{asin}$ as being *nearly* hopeless functions[^approx] and given the range of inputs we don't need to yet distinguish between one and two parameter $\text{atan}$.  There are a fair number of ways we can express $\theta$ via $\text{atan}$ and here are three of them:

$$
\begin{eqnarray}
\theta & = & \text{atan}\left( \frac{\sqrt{1-d^2}}{d} \right) \label{acos} \\
       & = & \frac{\pi}{2} - \text{ atan}\left( \frac{d}{\sqrt{1-d^2}} \right) \label{asin} \\
       & = & 2 \text{ atan}\left( \frac{\sqrt{1-d^2}}{1+d} \right) \label{hacos} 
\end{eqnarray}
$$

\\
where $\eqref{acos}$ and $\eqref{asin}$ are simply $\text{acos}$ and $\text{asin}$ in $\text{atan}$ form and $\eqref{hacos}$ is half-angle applied to $\eqref{acos}$.

Now for some practical considerations. By *reference* we can go two ways: close enough to detect problems in actually used methods or minimal error with respect to representable so we can accurately measure errors.  I'm going with the former.  Also we have to decide the meaning of a negative dot product.

1. It is an error so assert or spew a message.
2. Can I haz minimal anyway plz?
3. Go the long way around.

I like memes as much as the next person so will go with option 2.

Using C standard functions we could make a concise[^atan2sgn] implemention:

{% highlight c %}
void slerp_ref_0(quat_t* R, quat_t* A, quat_t* B, float t)
{
  float d   = quat_dot(A,B);    // cosine of angle between A & B = scalar part of BA^*
  float sgn = d >= 0 ? 1 : -1;  // need sign of 'd' for choosing minimal path
  float s0,s1;

  d = fabs(d);
  
  // deals with degenerate case
  if (d < SLERP_CUT) {
    // slerp part
    float s = sqrtf(1-d*d);    // sine of relative angle
    float a = atan2f(s,d);     // acos->atan2. some number of divides here :(
    float c = cosf(t*a);       // range reduction :(
    s1 = sqrtf(1-c*c)/s;       // atan2 might be computing 1/s :(
    s0 = c-d*s1;
  } else {
    // lerp
    s0 = 1.0f - t; 
    s1 = t;
  }

  // weighted sum: R = s0*A + sgn*s1*B
  quat_wsum(R, A, B, s0, sgn*s1);
}
{% endhighlight %}

\\
This leaves a fair number of thing we can know to the compiler to deduce and hopefully do something about.  Let's convert the *"slerp part"* to single parameter $\text{atan}$ with input on $\left[ 0,1 \right]$ using $\eqref{hacos}$:

{% highlight c %}
    // toy code: slerp_ref_1
    t += t;                      // account for half-angle atan
    float s2 = 1.f-d*d;
    float i  = recip_nr(1.f+d);  // ~x^-1  to at least 12-bit, 1 NR step
    float rs = rsqrt_nr(s2);     // ~x^-.5 to at least 12-bit, 1 NR step
    float y  = s2*rs;          
    float a  = atanf(y*i);       // still range reduction :(
    float s  = sinf(t*a);        // still range reduction :(
    float c  = sqrtf(1.f-s*s);   // can't ~1/sqrt(x) if max angle allowed
    s1 = s*rs;                   // sgn*s1*B can be computed while c in progress
    s0 = c-d*s1;
{% endhighlight %}

\\
Showing the use of approximation plus one step of newton-raphson for $x^{-\frac{1}{2}}$ and $x^{-1}$ since it has minimal impact on accuracy. Good idea or not is a different story.  Likewise for removing either of the forward trig ops.

Let's inject some more polynomial approximations:

{% highlight c %}
    // toy code: slerp_ref_3

    // set-up for atan
    float s2 = 1.f-d*d;
    float i  = recip_nr(1.f+d);
    float rs = rsqrt_nr(s2);
    float y  = s2*rs;
    float a  = atan_p(y*i);       // ~atan x on [0,1]
    
    // forward trig approx (or root for one)
    float ta = t*a;
    float c  = cos_4_5(ta);       // ~cos on [0, pi/4]
    float s  = sin_4_5(ta);       // ~sin on [0, pi/4]

    // double angle
    float tc = c+c;
    float cx = tc*c-1.f;          // double angle  (c^2-s^2 = 2c^2-1 = 1-2s^2)
    float sx = tc*s;              // double angle

    s1 = sx*rs;
    s0 = cx-d*s1;
{% endhighlight %}

\\
This changes to performing double-angle after the forward trig function(s) which allows less terms in these approximations.  Both are needed before we can continue so I'm showing performing both.  I'm not going to bother with trade-off comments anymore and will assume it's understood that I'm attempting to show various possible choices.

The only difficult to approximate function (with a small error) left is the arctangent. From the sum-of-tangents identity:

$$
\text{atan}\left(a\right) + \text{atan}\left(b\right)
= \text{atan}\left(\frac{a+b}{1-ab}\right)
$$

\\
To halve our angle range set $b=1$ since $\text{atan}\left(1\right)=\frac{\pi}{4}$.  Giving us the following:

$$
\frac{1}{2}\theta = \begin{cases}
\text{ atan}\left( \frac{\sqrt{1-d^2}-1}{d}\right)+\frac{\pi}{4}  & d \in
\left[0, \frac{\sqrt{2}}{2}\right) \\
\text{ atan}\left( \frac{\sqrt{1-d^2}}{1+d} \right) & d \in
\left[ \frac{\sqrt{2}}{2}, 1 \right] \\

\end{cases}
$$

\\
which limits the input to arctangent to $\left[0, \sqrt{2}-1\right]$

{% highlight c %}
  // toy code: slerp_ref_4
  float x  = (d < SQRT2_O_2) ? 1 : 0;  // select the range
  float ps = x*(PI*.25);               // multiply or select pi/4
  float s2 = 1.f-d*d;
  float i  = recip_nr(1.f-x+d);        // 1/d or 1/(1+d)
  float rs = rsqrt_nr(s2);   
  float y  = s2*rs-x;                  // sqrt(1-d^2)-1 or sqrt(1-d^2)
  float a  = atan_h_6(y*i)+ps;
{% endhighlight %}

\\
Sadly this last version explodes when 'd' is approaching zero.  This is easy enough to fix but this section is already getting a bit ridiculous..and for that matter 'd' approaching zero is pretty ridiculous as well.  Let's move on.

<br>

------

Precomputation <small>a small tweak to reference</small>
------

\\
A significant portion of the computational complexity in the previous section is computing values which are constant for fixed end-points.  Going back to $\eqref{slerpxy}$ we could instead precompute[^acos]:

$$ \begin{align*}
\theta & = \text{acos}\left(\left|d\right|\right) \\
X & = A \\
Y & = \text{sgn}\left(d\right)\frac{B - dA}{\sqrt{1 -d^2}}
\end{align*} $$

\\
which leaves using some method to calculate/approximate the two trig operations as before.  All the various computation choices from the reference section apply here as well so I'll just show the short version using standard functions:

{% highlight c %}
typedef struct {
  quat_t x,y;    // the 2D orthonormal basis
  float  a;      // the angle (+1 float to data-set)
} slerp_pc_t;

void slerp_pc_init_0(slerp_pc_t* k, quat_t* A, quat_t* B)
{
  float d   = quat_dot(A,B);
  float sgn = d >= 0 ? 1 : -1;

  d = fabs(d);
  quat_dup(&k->x, A);

  if (d < SLERP_CUT) {
    float s  = sqrtf(1-d*d);
    float a  = atan2f(s,d);
    float is = sgn/s;
    quat_wsum(&k->y, A, B, -d*is, is);
    k->a = b;
  }
  else {
    // can't be bothered to do this the right way
    quat_dup(&k->y, B);
    k->a = PI*0.5f;
  }
}

void slerp_pc_0(quat_t* r, slerp_pc_t* k, float t)
{
  float ta = t*k->a;

  // all of the previous options apply
  float c  = cosf(ta);
  float s  = sqrt(1-c*c);

  // or double the angle here

  // compute the weighted sum
  quat_wsum(r, &k->x, &k->y, c, s);
}
{% endhighlight %}

\\
Another option is replace using $A$ as the reference direction to the bisector of the end points:

$$ \begin{align*}
\theta & = \frac{1}{2} \text{acos}\left(\left|d\right|\right) \\
 X & = \left\Vert \text{sgn} (d) A + B\right\Vert \\
 Y & = \left\Vert B - \left(X \cdot B \right)X \right\Vert \\
\end{align*} $$

\\
and re-parameterize the input to $t' = t-.5$.  This transforms the range of $t\theta \in \left[-\frac{\pi}{4}, \frac{\pi}{4} \right]$, again to be polynomial approximation friendly.  This is equivalent to rotating the basis of the first option by $\frac{\theta}{2}$



<br>

------

Basic analysis tools <small></small>
------

\\
The action of slerp is defined by the principle real power $\eqref{power}$ and is bound to a complex plane.  So we can reason about it as walking a path[^notquite] from $\left(1,0\right) $ to $\left(a,b\right) = \left(\cos \theta,\sin \theta\right)$.  If the coordinate in the plane is described by:

$$ \begin{equation} \label{eq:fpoint}
p\left(t\right) = \left(p_x\!\left(t\right),~p_y\!\left(t\right)\right)
\end{equation} $$

\\
Then the angle at $t$ is:

$$ \begin{equation} \label{eq:fangle}
\phi\left(t\right) =\text{atan}\left(p_y\!\left(t\right), ~p_x\!\left(t\right)\right)
\end{equation} $$

\\
From $\phi\left(t\right)$ we can compute other quantites such as angular velocity $\left(\omega=\frac{d\phi}{dt}\right)$ and acceleration $\left(\alpha=\frac{d\omega}{dt}\right)$.

Taking the values for slerp $\eqref{power}$ and pluging them in we get the expected:

$$
\begin{eqnarray*}
  p\left(t\right) & = & \left( \cos\left(t\theta\right), ~\sin\left(t\theta\right)\right)  \\
  \phi\left(t\right) & = & \text{atan}\left( \sin\left(t\theta\right) , ~\cos\left(t\theta\right)\right) = t \theta \\
  \omega\left(t\right) & = & \theta \\
  \alpha\left(t\right) & = & 0
\end{eqnarray*}
$$

\\
This give a simple framework for computing errors which ignore the impact of floating-point computation.  We can also stick with a two-dimensional model for performing empirical testing

\\
Given some approximation of slerp we can compute the function $T$ to correct to constant angular velocity:

$$
\begin{eqnarray}
  \text{atan}\left(\frac{p_y\left(T(t)\right)}{p_x\left(T(t)\right)}\right) = t \theta \nonumber \\
  \frac{p_y\left(T(t)\right)}{p_x\left(T(t)\right)} = \tan\left(t \theta\right) \label{reparam} \\
\end{eqnarray}
$$

\\
by solving for $T$.

<br>

------

How good is LERP? <small> </small>
------

\\
If we approximate slerp with lerp then we have:

$$
\begin{eqnarray*}
  p_l\left(t\right)      & = & \left( 1+t\left(\cos\left(\theta\right)-1\right), ~t\sin\left(\theta\right)\right) \\
  \phi_l\left(t\right)   & = & \text{atan}\left(~t\sin\left(\theta\right), ~1+t\left(\cos\left(\theta\right)-1\right)\right) \\
  \omega_l\left(t\right) & = & \frac{\sin\left(\theta\right)}{1-2t\left(\cos\left(\theta\right) \left(t-1\right)  - \left(t-1\right)  \right)} \\
  \alpha_l\left(t\right) & = & \frac{2\sin\left(\theta\right)\left(\cos\left(\theta\right)-1\right)\left(2t-1\right)}
  {\left(1+2t\left(\cos\left(\theta\right) \left(t-1\right)  - \left(t-1\right)  \right)\right)^2}
\end{eqnarray*}
$$

ff

<div id="fig1" style="width:100%"></div>

ff

<div id="fig2" style="width:100%"></div>

<br>

------

Reparm LERP <small></small>
------

Since the b

$$
  \frac{r(t)~\sin \left( \theta \right) }{1+r(t) \left( \cos \left( \theta \right)-1 \right)}  =  \tan \left(t \theta \right) 
$$

\\
Solving for $r(t)$:


$$ \frac{1}{\sqrt{1-d^2} \cot \left(t~\text{acos}(d)\right)-d+1} $$

J Blow [^jblow] xxx Arseny Kapoulkine[^zeux]


<br>

------

Half-angle transformed input
------


<br>

------

Eberly's approximation
------

\\
David Eberly in 2011[^eberly]


<br>

------

Cayley transformed input
------


\\
asdfasdf

<br>

------

References and Footnotes
------

[^shoemake85]: *"Animating Rotation with Quaternion Curves"*, Ken Shoemake, 1985. ([PDF](http://run.usc.edu/cs520-s12/assign2/p245-shoemake.pdf))

[^notquite]:   Actually any scalar multiple of the end-points is mathematically correct.
[^nosign]:     no sign fix-up needed.
[^zeux]:       *"Approximating slerp"*, Arseny Kapoulkine, 2015 ([PAGE](http://zeuxcg.org/2015/07/23/approximating-slerp/))
[^jblow]:      *"Hacking Quaternions"*, Jonathan Blow, 2002 ([PAGE](http://number-none.com/product/Hacking%20Quaternions/))
[^approx]:     but might be worth consideration if building approximate functions
[^opt]:        I'm of the school of thought that its sinful to leave easy things to know to be true to the compiler to hopefully deduce and do something about.
[^eberly]:     *"A Fast and Accurate Algorithm for Computing SLERP"*, David Eberly, 2011.
[^id]:        *"Slerping Clock Cycles"*, van Waveren, 2005.
[^acos]:      Showing arccosine to be short, previous methods still apply
[^sincr]:     arctanget is strictly increasing on the interval
[^atan2sgn]:  The abs and sgn computations aren't needed with atan2.  Assuming nobody's really going to use the versions.


<script>
var plot0 = {
x:[2.04082e-8, 0.000306718, 0.000613415, 0.00122681, 0.0024536, 
0.00490718, 0.00981434, 0.0196287, 0.0409084, 0.0607779, 0.0802576, 
0.101388, 0.121109, 0.121443, 0.121777, 0.122445, 0.123781, 0.126452, 
0.131795, 0.142481, 0.142809, 0.143137, 0.143792, 0.145104, 0.147726, 
0.152972, 0.1533, 0.153627, 0.154283, 0.155595, 0.158217, 0.163463, 
0.163769, 0.164074, 0.164686, 0.165909, 0.168356, 0.173249, 0.173554, 
0.17386, 0.174472, 0.175695, 0.178142, 0.183035, 0.183366, 0.183698, 
0.184361, 0.185687, 0.18834, 0.188672, 0.189003, 0.189667, 0.190993, 
0.193646, 0.193978, 0.194309, 0.194972, 0.196299, 0.198952, 0.204257, 
0.204567, 0.204877, 0.205496, 0.206734, 0.209211, 0.20952, 0.20983, 
0.210449, 0.211687, 0.214164, 0.214473, 0.214783, 0.215402, 0.21664, 
0.21695, 0.217259, 0.217879, 0.219117, 0.219426, 0.219736, 0.220355, 
0.221593, 0.22407, 0.224374, 0.224677, 0.225284, 0.226498, 0.226801, 
0.227105, 0.227712, 0.228926, 0.229229, 0.229533, 0.23014, 0.230443, 
0.230747, 0.231354, 0.231657, 0.231961, 0.232568, 0.232871, 0.233175, 
0.233781, 0.234085, 0.234388, 0.234692, 0.234995, 0.235299, 0.235602, 
0.235906, 0.236209, 0.236513, 0.236816, 0.23712, 0.237423, 0.237727, 
0.23803, 0.238334, 0.238637, 0.238941, 0.239244, 0.239548, 0.239851, 
0.240155, 0.240458, 0.240762, 0.241065, 0.241369, 0.241672, 0.241976, 
0.242279, 0.242582, 0.242886, 0.243493, 0.243822, 0.244151, 0.24481, 
0.245139, 0.245469, 0.246127, 0.246456, 0.246786, 0.247444, 0.248761, 
0.249091, 0.24942, 0.250079, 0.251396, 0.251725, 0.252054, 0.252713, 
0.25403, 0.254359, 0.254688, 0.255347, 0.256664, 0.256993, 0.257323, 
0.257981, 0.259298, 0.259628, 0.259957, 0.260616, 0.261933, 0.264567, 
0.264874, 0.265181, 0.265796, 0.267025, 0.269483, 0.26979, 0.270097, 
0.270712, 0.271941, 0.274399, 0.274706, 0.275013, 0.275628, 0.276857, 
0.279315, 0.284231, 0.284564, 0.284897, 0.285563, 0.286895, 0.289559, 
0.294888, 0.295221, 0.295554, 0.29622, 0.297552, 0.300217, 0.305546, 
0.305872, 0.306199, 0.306853, 0.308161, 0.310777, 0.316008, 0.326471, 
0.326776, 0.32708, 0.32769, 0.32891, 0.331349, 0.336228, 0.345985, 
0.367151, 0.386907, 0.408314, 0.429331, 0.448938, 0.470196, 0.490044, 
0.509502, 0.530611, 0.55031, 0.57166, 0.59262, 0.61217, 0.633371, 
0.653162, 0.653465, 0.653769, 0.654375, 0.655587, 0.658013, 0.662863, 
0.672563, 0.672892, 0.673221, 0.673879, 0.675195, 0.677827, 0.68309, 
0.683419, 0.683747, 0.684405, 0.685721, 0.688353, 0.693616, 0.693923, 
0.694229, 0.694843, 0.696071, 0.698526, 0.703437, 0.703744, 0.704051, 
0.704664, 0.705892, 0.708347, 0.713258, 0.71359, 0.713923, 0.714589, 
0.715919, 0.718581, 0.718914, 0.719246, 0.719912, 0.721243, 0.723904, 
0.724237, 0.72457, 0.725235, 0.726566, 0.729228, 0.72956, 0.729893, 
0.730558, 0.731889, 0.734551, 0.734862, 0.735172, 0.735794, 0.737036, 
0.739522, 0.739832, 0.740143, 0.740764, 0.742007, 0.744492, 0.744803, 
0.745114, 0.745735, 0.746978, 0.747288, 0.747599, 0.74822, 0.749463, 
0.749774, 0.750084, 0.750706, 0.751016, 0.751327, 0.751948, 0.752259, 
0.75257, 0.753191, 0.753502, 0.753812, 0.754434, 0.754738, 0.755043, 
0.755652, 0.755957, 0.756261, 0.75687, 0.757175, 0.75748, 0.757784, 
0.758089, 0.758393, 0.758698, 0.759003, 0.759307, 0.759612, 0.759916, 
0.760221, 0.760525, 0.76083, 0.761135, 0.761439, 0.761744, 0.762048, 
0.762353, 0.762658, 0.762962, 0.763267, 0.763571, 0.763876, 0.76418, 
0.764485, 0.76479, 0.765094, 0.765399, 0.765703, 0.766008, 0.766617, 
0.766922, 0.767226, 0.767835, 0.76814, 0.768445, 0.769054, 0.769358, 
0.769663, 0.770272, 0.77149, 0.771795, 0.7721, 0.772709, 0.773927, 
0.774257, 0.774588, 0.775249, 0.77657, 0.7769, 0.777231, 0.777892, 
0.779213, 0.779543, 0.779874, 0.780535, 0.781856, 0.784499, 0.78483, 
0.78516, 0.785821, 0.787142, 0.789785, 0.790116, 0.790446, 0.791107, 
0.792428, 0.795071, 0.79538, 0.795688, 0.796305, 0.797538, 0.800005, 
0.804938, 0.805247, 0.805555, 0.806172, 0.807405, 0.809872, 0.814805, 
0.815139, 0.815474, 0.816142, 0.817478, 0.820152, 0.825498, 0.825832, 
0.826166, 0.826834, 0.828171, 0.830844, 0.83619, 0.836519, 0.836847, 
0.837503, 0.838815, 0.841439, 0.846688, 0.857186, 0.857492, 0.857798, 
0.85841, 0.859634, 0.862082, 0.866978, 0.876771, 0.898007, 0.917833, 
0.93727, 0.958357, 0.978034, 0.978377, 0.978721, 0.979407, 0.98078, 
0.983526, 0.989017, 0.98936, 0.989704, 0.99039, 0.991763, 0.994509, 
0.994852, 0.995195, 0.995881, 0.997254, 0.997597, 0.997941, 0.998627, 
0.99897, 0.999314, 0.999657, 1.],
y:[3.70796e-9, 0.0000556977, 0.000111332, 0.00022242, 0.000443876, 
0.000883896, 0.00175231, 0.00344209, 0.00688546, 0.00981951, 
0.012423, 0.0149312, 0.0169667, 0.0169986, 0.0170304, 0.0170937, 
0.0172193, 0.0174663, 0.0179436, 0.0188305, 0.0188563, 0.018882, 
0.0189331, 0.0190343, 0.0192327, 0.0196127, 0.0196357, 0.0196587, 
0.0197043, 0.0197945, 0.0199707, 0.0203063, 0.0203252, 0.020344, 
0.0203814, 0.0204553, 0.0205993, 0.0208728, 0.0208892, 0.0209056, 
0.0209381, 0.0210022, 0.0211267, 0.021361, 0.0213762, 0.0213913, 
0.0214211, 0.0214798, 0.0215929, 0.0216066, 0.0216202, 0.0216472, 
0.0217001, 0.0218015, 0.0218138, 0.021826, 0.0218501, 0.0218972, 
0.021987, 0.0221492, 0.0221579, 0.0221666, 0.0221837, 0.0222169, 
0.0222796, 0.0222871, 0.0222945, 0.022309, 0.0223372, 0.0223897, 
0.022396, 0.0224021, 0.0224141, 0.0224372, 0.0224428, 0.0224483, 
0.022459, 0.0224796, 0.0224845, 0.0224894, 0.0224989, 0.0225169, 
0.0225492, 0.0225528, 0.0225563, 0.0225631, 0.0225759, 0.0225789, 
0.0225818, 0.0225874, 0.0225977, 0.0226001, 0.0226024, 0.0226068, 
0.0226089, 0.0226109, 0.0226147, 0.0226165, 0.0226182, 0.0226213, 
0.0226228, 0.0226242, 0.0226268, 0.022628, 0.0226291, 0.0226301, 
0.022631, 0.0226319, 0.0226327, 0.0226334, 0.0226341, 0.0226346, 
0.0226351, 0.0226355, 0.0226359, 0.0226362, 0.0226363, 0.0226365, 
0.0226365, 0.0226365, 0.0226363, 0.0226362, 0.0226359, 0.0226355, 
0.0226351, 0.0226346, 0.0226341, 0.0226334, 0.0226327, 0.0226319, 
0.0226311, 0.0226301, 0.0226291, 0.0226268, 0.0226255, 0.022624, 
0.0226209, 0.0226192, 0.0226174, 0.0226135, 0.0226114, 0.0226093, 
0.0226047, 0.0225945, 0.0225917, 0.0225889, 0.0225829, 0.0225699, 
0.0225664, 0.0225628, 0.0225554, 0.0225395, 0.0225354, 0.0225311, 
0.0225223, 0.0225036, 0.0224987, 0.0224937, 0.0224835, 0.022462, 
0.0224565, 0.0224508, 0.0224392, 0.0224149, 0.0223621, 0.0223556, 
0.022349, 0.0223356, 0.0223079, 0.0222488, 0.022241, 0.0222332, 
0.0222174, 0.0221849, 0.0221161, 0.0221072, 0.0220982, 0.02208, 
0.0220426, 0.0219643, 0.0217935, 0.0217812, 0.0217689, 0.0217439, 
0.021693, 0.021587, 0.0213585, 0.0213435, 0.0213284, 0.021298, 
0.0212361, 0.0211083, 0.0208367, 0.0208193, 0.0208019, 0.0207667, 
0.0206955, 0.0205493, 0.0202418, 0.0195674, 0.0195466, 0.0195257, 
0.0194837, 0.019399, 0.0192264, 0.018869, 0.0181058, 0.0162398, 
0.0142593, 0.0118848, 0.00935872, 0.00686376, 0.00404938, 0.00135936, 
-0.00129742, -0.00415778, -0.0067659, -0.0094816, -0.0119929, 
-0.0141616, -0.0162892, -0.0180362, -0.018061, -0.0180858, 
-0.0181351, -0.0182331, -0.0184262, -0.0188006, -0.0195012, 
-0.0195238, -0.0195464, -0.0195912, -0.01968, -0.0198539, -0.0201868, 
-0.0202069, -0.020227, -0.0202668, -0.0203457, -0.0204995, -0.020792, 
-0.0208084, -0.0208247, -0.0208572, -0.0209213, -0.0210462, 
-0.0212822, -0.0212963, -0.0213104, -0.0213383, -0.0213933, 
-0.0214998, -0.0216989, -0.0217117, -0.0217244, -0.0217496, 
-0.021799, -0.0218935, -0.0219049, -0.0219163, -0.0219387, 
-0.0219825, -0.0220659, -0.0220759, -0.0220859, -0.0221055, 
-0.0221437, -0.0222158, -0.0222245, -0.022233, -0.0222498, 
-0.0222824, -0.0223432, -0.0223499, -0.0223566, -0.0223697, 
-0.0223949, -0.0224416, -0.0224471, -0.0224525, -0.0224631, 
-0.0224833, -0.0225201, -0.0225243, -0.0225285, -0.0225366, 
-0.0225518, -0.0225554, -0.022559, -0.0225658, -0.0225785, 
-0.0225815, -0.0225844, -0.02259, -0.0225926, -0.0225952, -0.0226002, 
-0.0226025, -0.0226048, -0.0226091, -0.0226112, -0.0226131, 
-0.0226168, -0.0226185, -0.0226201, -0.0226231, -0.0226245, 
-0.0226258, -0.0226282, -0.0226293, -0.0226303, -0.0226312, 
-0.0226321, -0.0226329, -0.0226336, -0.0226342, -0.0226348, 
-0.0226352, -0.0226356, -0.022636, -0.0226362, -0.0226364, 
-0.0226365, -0.0226365, -0.0226364, -0.0226363, -0.0226361, 
-0.0226358, -0.0226354, -0.022635, -0.0226345, -0.0226339, 
-0.0226332, -0.0226325, -0.0226317, -0.0226308, -0.0226298, 
-0.0226287, -0.0226276, -0.0226251, -0.0226238, -0.0226224, 
-0.0226193, -0.0226176, -0.0226159, -0.0226122, -0.0226102, 
-0.0226082, -0.0226038, -0.0225943, -0.0225917, -0.0225891, 
-0.0225835, -0.0225715, -0.0225681, -0.0225645, -0.0225571, 
-0.0225413, -0.0225371, -0.0225329, -0.022524, -0.0225053, 
-0.0225004, -0.0224954, -0.0224851, -0.0224635, -0.022416, 
-0.0224096, -0.0224032, -0.0223901, -0.0223627, -0.0223036, 
-0.0222958, -0.0222879, -0.0222719, -0.0222387, -0.022168, 
-0.0221594, -0.0221507, -0.0221331, -0.0220969, -0.0220207, 
-0.0218533, -0.0218422, -0.0218309, -0.0218083, -0.021762, 
-0.0216658, -0.0214582, -0.0214434, -0.0214285, -0.0213985, 
-0.0213374, -0.0212107, -0.0209397, -0.020922, -0.0209042, 
-0.0208683, -0.0207955, -0.0206454, -0.0203277, -0.0203075, 
-0.0202872, -0.0202462, -0.0201633, -0.0199933, -0.0196366, 
-0.0188567, -0.0188327, -0.0188086, -0.0187601, -0.0186622, 
-0.0184629, -0.01805, -0.0171676, -0.014998, -0.0126632, -0.0100928, 
-0.00699883, -0.00383513, -0.00377763, -0.00372006, -0.00360469, 
-0.00337299, -0.00290588, -0.0019568, -0.00189683, -0.00183678, 
-0.00171646, -0.00147489, -0.000988108, -0.000926918, -0.000865652, 
-0.000742893, -0.000496469, -0.000434674, -0.000372804, -0.000248838, 
-0.000186742, -0.000124571, -0.0000623248, -3.70796e-9],
name: 'pi/2',
type: 'scatter'
};

var plot1 = {
x:[2.04082e-8, 0.000306718, 0.000613415, 0.00122681, 0.0024536, 
0.00490718, 0.00981434, 0.0196287, 0.0409084, 0.0607779, 0.0802576, 
0.0805878, 0.080918, 0.0815783, 0.082899, 0.0855404, 0.0908231, 
0.101388, 0.101697, 0.102005, 0.102621, 0.103854, 0.106319, 0.111249, 
0.111557, 0.111865, 0.112481, 0.113714, 0.116179, 0.121109, 0.121443, 
0.121777, 0.122445, 0.123781, 0.126452, 0.131795, 0.132129, 0.132463, 
0.133131, 0.134466, 0.137138, 0.142481, 0.142809, 0.143137, 0.143792, 
0.145104, 0.147726, 0.152972, 0.1533, 0.153627, 0.154283, 0.155595, 
0.158217, 0.163463, 0.163769, 0.164074, 0.164686, 0.165909, 0.168356, 
0.168662, 0.168967, 0.169579, 0.170802, 0.173249, 0.173554, 0.17386, 
0.174472, 0.175695, 0.178142, 0.183035, 0.183366, 0.183698, 0.184361, 
0.185687, 0.18834, 0.188672, 0.189003, 0.189667, 0.190993, 0.193646, 
0.193978, 0.194309, 0.194972, 0.196299, 0.19663, 0.196962, 0.197625, 
0.198952, 0.199283, 0.199615, 0.200278, 0.201605, 0.201936, 0.202268, 
0.202931, 0.204257, 0.204567, 0.204877, 0.205496, 0.205805, 0.206115, 
0.206734, 0.207044, 0.207353, 0.207972, 0.209211, 0.20952, 0.20983, 
0.210449, 0.210758, 0.211068, 0.211687, 0.211997, 0.212306, 0.212616, 
0.212925, 0.213235, 0.213545, 0.213854, 0.214164, 0.214473, 0.214783, 
0.215092, 0.215402, 0.215712, 0.216021, 0.216331, 0.21664, 0.21695, 
0.217259, 0.217569, 0.217879, 0.218188, 0.218498, 0.218807, 0.219117, 
0.219426, 0.219736, 0.220046, 0.220355, 0.220665, 0.220974, 0.221284, 
0.221593, 0.221903, 0.222213, 0.222832, 0.223141, 0.223451, 0.22407, 
0.224374, 0.224677, 0.225284, 0.225587, 0.225891, 0.226498, 0.226801, 
0.227105, 0.227712, 0.228926, 0.229229, 0.229533, 0.23014, 0.231354, 
0.233781, 0.234085, 0.234388, 0.234995, 0.236209, 0.238637, 0.238941, 
0.239244, 0.239851, 0.241065, 0.243493, 0.243822, 0.244151, 0.24481, 
0.246127, 0.248761, 0.249091, 0.24942, 0.250079, 0.251396, 0.25403, 
0.254359, 0.254688, 0.255347, 0.256664, 0.259298, 0.264567, 0.264874, 
0.265181, 0.265796, 0.267025, 0.269483, 0.274399, 0.274706, 0.275013, 
0.275628, 0.276857, 0.279315, 0.284231, 0.284564, 0.284897, 0.285563, 
0.286895, 0.289559, 0.294888, 0.305546, 0.305872, 0.306199, 0.306853, 
0.308161, 0.310777, 0.316008, 0.326471, 0.345985, 0.367151, 0.386907, 
0.408314, 0.429331, 0.448938, 0.470196, 0.490044, 0.509502, 0.530611, 
0.55031, 0.57166, 0.59262, 0.61217, 0.633371, 0.653162, 0.672563, 
0.672892, 0.673221, 0.673879, 0.675195, 0.677827, 0.68309, 0.693616, 
0.693923, 0.694229, 0.694843, 0.696071, 0.698526, 0.703437, 0.703744, 
0.704051, 0.704664, 0.705892, 0.708347, 0.713258, 0.71359, 0.713923, 
0.714589, 0.715919, 0.718581, 0.723904, 0.724237, 0.72457, 0.725235, 
0.726566, 0.729228, 0.734551, 0.734862, 0.735172, 0.735794, 0.737036, 
0.739522, 0.739832, 0.740143, 0.740764, 0.742007, 0.744492, 0.744803, 
0.745114, 0.745735, 0.746978, 0.749463, 0.754434, 0.754738, 0.755043, 
0.755652, 0.75687, 0.759307, 0.759612, 0.759916, 0.760525, 0.761744, 
0.76418, 0.764485, 0.76479, 0.765399, 0.766617, 0.766922, 0.767226, 
0.767835, 0.769054, 0.769358, 0.769663, 0.770272, 0.77149, 0.771795, 
0.7721, 0.772709, 0.773927, 0.774257, 0.774588, 0.775249, 0.775579, 
0.775909, 0.77657, 0.7769, 0.777231, 0.777561, 0.777892, 0.778222, 
0.778552, 0.778883, 0.779213, 0.779543, 0.779874, 0.780204, 0.780535, 
0.780865, 0.781195, 0.781526, 0.781856, 0.782186, 0.782517, 0.782847, 
0.783178, 0.783508, 0.783838, 0.784169, 0.784499, 0.78483, 0.78516, 
0.78549, 0.785821, 0.786151, 0.786481, 0.787142, 0.787473, 0.787803, 
0.788464, 0.788794, 0.789124, 0.789785, 0.790116, 0.790446, 0.791107, 
0.791437, 0.791767, 0.792428, 0.792759, 0.793089, 0.79375, 0.795071, 
0.79538, 0.795688, 0.796305, 0.797538, 0.797846, 0.798155, 0.798771, 
0.800005, 0.800313, 0.800621, 0.801238, 0.802472, 0.804938, 0.805247, 
0.805555, 0.806172, 0.807405, 0.809872, 0.81018, 0.810488, 0.811105, 
0.812339, 0.814805, 0.815139, 0.815474, 0.816142, 0.817478, 0.820152, 
0.820486, 0.82082, 0.821488, 0.822825, 0.825498, 0.825832, 0.826166, 
0.826834, 0.828171, 0.830844, 0.83619, 0.836519, 0.836847, 0.837503, 
0.838815, 0.841439, 0.846688, 0.847016, 0.847344, 0.848, 0.849313, 
0.851937, 0.857186, 0.857492, 0.857798, 0.85841, 0.859634, 0.862082, 
0.866978, 0.876771, 0.877103, 0.877435, 0.878098, 0.879426, 0.88208, 
0.887389, 0.898007, 0.898317, 0.898627, 0.899246, 0.900486, 0.902964, 
0.90792, 0.917833, 0.93727, 0.958357, 0.978034, 0.978377, 0.978721, 
0.979407, 0.98078, 0.983526, 0.989017, 0.98936, 0.989704, 0.99039, 
0.991763, 0.994509, 0.994852, 0.995195, 0.995881, 0.997254, 0.997597, 
0.997941, 0.998627, 0.99897, 0.999314, 0.999657, 1.],
y:[5.0859e-10, 7.63749e-6, 0.0000152621, 0.0000304741, 0.0000607494, 
0.000120706, 0.00023825, 0.000463907, 0.000910474, 0.00127558, 
0.00158602, 0.00159088, 0.00159573, 0.00160539, 0.00162456, 
0.00166225, 0.00173513, 0.00187094, 0.0018747, 0.00187845, 
0.00188592, 0.00190073, 0.00192982, 0.00198587, 0.00198928, 
0.00199268, 0.00199944, 0.00201284, 0.00203911, 0.00208957, 
0.00209288, 0.00209619, 0.00210276, 0.00211575, 0.00214113, 
0.00218948, 0.0021924, 0.0021953, 0.00220107, 0.00221246, 0.00223464, 
0.00227666, 0.00227913, 0.0022816, 0.00228649, 0.00229613, 
0.00231486, 0.00235009, 0.0023522, 0.00235429, 0.00235844, 0.0023666, 
0.00238237, 0.00241174, 0.00241336, 0.00241497, 0.00241817, 
0.00242444, 0.00243652, 0.00243799, 0.00243945, 0.00244233, 
0.00244798, 0.00245883, 0.00246014, 0.00246144, 0.00246402, 
0.00246906, 0.00247867, 0.00249609, 0.00249718, 0.00249826, 
0.00250039, 0.00250452, 0.00251225, 0.00251317, 0.00251408, 
0.00251586, 0.00251929, 0.00252563, 0.00252638, 0.00252711, 
0.00252854, 0.00253128, 0.00253194, 0.00253259, 0.00253385, 
0.00253625, 0.00253683, 0.00253739, 0.00253848, 0.00254054, 
0.00254103, 0.00254151, 0.00254243, 0.00254416, 0.00254453, 
0.0025449, 0.00254561, 0.00254595, 0.00254629, 0.00254692, 
0.00254723, 0.00254753, 0.00254809, 0.00254911, 0.00254935, 
0.00254957, 0.00254999, 0.00255019, 0.00255038, 0.00255072, 
0.00255089, 0.00255104, 0.00255118, 0.00255132, 0.00255144, 
0.00255156, 0.00255166, 0.00255176, 0.00255185, 0.00255193, 
0.00255201, 0.00255207, 0.00255212, 0.00255217, 0.00255221, 
0.00255223, 0.00255225, 0.00255226, 0.00255226, 0.00255226, 
0.00255224, 0.00255221, 0.00255218, 0.00255214, 0.00255209, 
0.00255203, 0.00255196, 0.00255188, 0.00255179, 0.0025517, 
0.00255159, 0.00255148, 0.00255136, 0.00255123, 0.00255094, 
0.00255079, 0.00255062, 0.00255027, 0.00255008, 0.00254989, 
0.00254947, 0.00254925, 0.00254902, 0.00254854, 0.00254829, 
0.00254803, 0.00254748, 0.00254629, 0.00254597, 0.00254564, 
0.00254496, 0.00254351, 0.00254021, 0.00253976, 0.0025393, 
0.00253837, 0.00253639, 0.00253206, 0.00253149, 0.0025309, 
0.00252971, 0.00252723, 0.00252189, 0.00252112, 0.00252035, 
0.00251878, 0.00251552, 0.00250858, 0.00250767, 0.00250675, 
0.00250488, 0.00250105, 0.00249295, 0.00249189, 0.00249083, 
0.00248868, 0.00248428, 0.00247504, 0.0024549, 0.00245365, 0.0024524, 
0.00244988, 0.00244475, 0.00243412, 0.00241147, 0.00240999, 
0.00240851, 0.00240551, 0.00239945, 0.00238697, 0.00236067, 
0.00235882, 0.00235697, 0.00235323, 0.00234567, 0.00233016, 
0.00229761, 0.00222658, 0.00222428, 0.00222197, 0.00221734, 
0.00220798, 0.00218892, 0.00214947, 0.00206539, 0.00189121, 
0.00167915, 0.00146218, 0.00120945, 0.000946672, 0.000691349, 
0.000406614, 0.000136305, -0.000130092, -0.000417536, -0.000681399, 
-0.000959342, -0.00122081, -0.00145164, -0.00168465, -0.00188312, 
-0.00205729, -0.00206005, -0.00206281, -0.00206831, -0.00207922, 
-0.00210074, -0.00214248, -0.00222066, -0.00222284, -0.002225, 
-0.00222931, -0.00223785, -0.00225463, -0.00228696, -0.00228893, 
-0.00229089, -0.00229478, -0.0023025, -0.00231761, -0.00234654, 
-0.00234844, -0.00235033, -0.00235409, -0.0023615, -0.00237593, 
-0.00240322, -0.00240485, -0.00240648, -0.00240971, -0.00241606, 
-0.00242836, -0.00245131, -0.00245258, -0.00245384, -0.00245634, 
-0.00246126, -0.00247072, -0.00247187, -0.00247301, -0.00247527, 
-0.00247969, -0.00248815, -0.00248918, -0.00249019, -0.0024922, 
-0.00249611, -0.00250357, -0.00251693, -0.00251768, -0.00251842, 
-0.00251989, -0.00252272, -0.002528, -0.00252863, -0.00252924, 
-0.00253045, -0.00253278, -0.00253704, -0.00253754, -0.00253803, 
-0.00253898, -0.00254079, -0.00254122, -0.00254164, -0.00254246, 
-0.00254401, -0.00254438, -0.00254473, -0.00254543, -0.00254671, 
-0.00254701, -0.0025473, -0.00254786, -0.00254888, -0.00254914, 
-0.00254938, -0.00254984, -0.00255005, -0.00255026, -0.00255063, 
-0.00255081, -0.00255097, -0.00255113, -0.00255128, -0.00255141, 
-0.00255154, -0.00255165, -0.00255176, -0.00255185, -0.00255194, 
-0.00255201, -0.00255208, -0.00255214, -0.00255218, -0.00255222, 
-0.00255224, -0.00255226, -0.00255226, -0.00255226, -0.00255225, 
-0.00255222, -0.00255219, -0.00255214, -0.00255209, -0.00255202, 
-0.00255195, -0.00255186, -0.00255177, -0.00255166, -0.00255155, 
-0.00255129, -0.00255114, -0.00255099, -0.00255064, -0.00255046, 
-0.00255026, -0.00254984, -0.00254961, -0.00254937, -0.00254887, 
-0.0025486, -0.00254832, -0.00254773, -0.00254742, -0.0025471, 
-0.00254643, -0.00254496, -0.0025446, -0.00254422, -0.00254345, 
-0.00254178, -0.00254135, -0.0025409, -0.00253998, -0.00253802, 
-0.00253751, -0.00253699, -0.00253592, -0.00253367, -0.00252873, 
-0.00252807, -0.00252741, -0.00252604, -0.0025232, -0.00251707, 
-0.00251626, -0.00251545, -0.00251378, -0.00251034, -0.00250301, 
-0.00250197, -0.00250092, -0.00249878, -0.00249437, -0.00248502, 
-0.0024838, -0.00248257, -0.00248007, -0.00247495, -0.00246414, 
-0.00246274, -0.00246133, -0.00245847, -0.00245261, -0.00244034, 
-0.00241358, -0.00241184, -0.00241008, -0.00240655, -0.00239934, 
-0.00238438, -0.00235227, -0.00235017, -0.00234805, -0.00234379, 
-0.00233512, -0.00231721, -0.00227917, -0.00227686, -0.00227454, 
-0.00226987, -0.0022604, -0.00224097, -0.00220013, -0.00211041, 
-0.00210718, -0.00210394, -0.00209742, -0.00208423, -0.00205724, 
-0.00200086, -0.00187831, -0.00187453, -0.00187075, -0.00186315, 
-0.00184781, -0.00181659, -0.00175197, -0.00161396, -0.0013088, 
-0.000924857, -0.000515798, -0.000508221, -0.00050063, -0.000485401, 
-0.000454761, -0.000392746, -0.000265776, -0.00025771, -0.000249628, 
-0.000233419, -0.000200815, -0.000134867, -0.000126554, -0.000118225, 
-0.000101522, -0.0000679296, -0.0000594927, -0.0000510404, 
-0.0000340892, -0.0000255903, -0.0000170759, -8.54596e-6, -5.0859e-10],
mode: 'lines',
name: 'pi/4'
};

var plot2 = {
x:[2.04082e-8, 0.000306718, 0.000613415, 0.00122681, 0.0024536, 
0.00490718, 0.00981434, 0.0196287, 0.0409084, 0.0607779, 0.0610823, 
0.0613866, 0.0619954, 0.0632129, 0.0656478, 0.0705178, 0.0802576, 
0.0805878, 0.080918, 0.0815783, 0.082899, 0.0855404, 0.0908231, 
0.101388, 0.101697, 0.102005, 0.102621, 0.103854, 0.106319, 0.111249, 
0.111557, 0.111865, 0.112481, 0.113714, 0.116179, 0.121109, 0.121443, 
0.121777, 0.122445, 0.123781, 0.126452, 0.131795, 0.132129, 0.132463, 
0.133131, 0.134466, 0.137138, 0.142481, 0.142809, 0.143137, 0.143792, 
0.145104, 0.147726, 0.152972, 0.1533, 0.153627, 0.154283, 0.155595, 
0.158217, 0.163463, 0.163769, 0.164074, 0.164686, 0.165909, 0.168356, 
0.168662, 0.168967, 0.169579, 0.170802, 0.173249, 0.173554, 0.17386, 
0.174472, 0.175695, 0.178142, 0.178447, 0.178753, 0.179365, 0.180588, 
0.183035, 0.183366, 0.183698, 0.184361, 0.185687, 0.18834, 0.188672, 
0.189003, 0.189667, 0.190993, 0.193646, 0.193978, 0.194309, 0.194972, 
0.196299, 0.19663, 0.196962, 0.197625, 0.198952, 0.199283, 0.199615, 
0.200278, 0.20061, 0.200941, 0.201605, 0.201936, 0.202268, 0.202931, 
0.204257, 0.204567, 0.204877, 0.205496, 0.205805, 0.206115, 0.206734, 
0.207044, 0.207353, 0.207663, 0.207972, 0.208282, 0.208591, 0.209211, 
0.20952, 0.20983, 0.210139, 0.210449, 0.210758, 0.211068, 0.211378, 
0.211687, 0.211997, 0.212306, 0.212616, 0.212925, 0.213235, 0.213545, 
0.213854, 0.214164, 0.214473, 0.214783, 0.215092, 0.215402, 0.215712, 
0.216021, 0.216331, 0.21664, 0.21695, 0.217259, 0.217569, 0.217879, 
0.218188, 0.218498, 0.219117, 0.219426, 0.219736, 0.220355, 0.220665, 
0.220974, 0.221593, 0.221903, 0.222213, 0.222832, 0.22407, 0.224374, 
0.224677, 0.225284, 0.226498, 0.226801, 0.227105, 0.227712, 0.228926, 
0.229229, 0.229533, 0.23014, 0.231354, 0.233781, 0.234085, 0.234388, 
0.234995, 0.236209, 0.238637, 0.238941, 0.239244, 0.239851, 0.241065, 
0.243493, 0.243822, 0.244151, 0.24481, 0.246127, 0.248761, 0.25403, 
0.254359, 0.254688, 0.255347, 0.256664, 0.259298, 0.264567, 0.264874, 
0.265181, 0.265796, 0.267025, 0.269483, 0.274399, 0.274706, 0.275013, 
0.275628, 0.276857, 0.279315, 0.284231, 0.284564, 0.284897, 0.285563, 
0.286895, 0.289559, 0.294888, 0.305546, 0.305872, 0.306199, 0.306853, 
0.308161, 0.310777, 0.316008, 0.326471, 0.345985, 0.367151, 0.386907, 
0.408314, 0.429331, 0.448938, 0.470196, 0.490044, 0.509502, 0.530611, 
0.55031, 0.57166, 0.59262, 0.61217, 0.633371, 0.653162, 0.672563, 
0.672892, 0.673221, 0.673879, 0.675195, 0.677827, 0.68309, 0.693616, 
0.693923, 0.694229, 0.694843, 0.696071, 0.698526, 0.703437, 0.703744, 
0.704051, 0.704664, 0.705892, 0.708347, 0.713258, 0.71359, 0.713923, 
0.714589, 0.715919, 0.718581, 0.723904, 0.724237, 0.72457, 0.725235, 
0.726566, 0.729228, 0.734551, 0.734862, 0.735172, 0.735794, 0.737036, 
0.739522, 0.744492, 0.744803, 0.745114, 0.745735, 0.746978, 0.749463, 
0.754434, 0.754738, 0.755043, 0.755652, 0.75687, 0.759307, 0.759612, 
0.759916, 0.760525, 0.761744, 0.76418, 0.764485, 0.76479, 0.765399, 
0.766617, 0.766922, 0.767226, 0.767835, 0.769054, 0.769358, 0.769663, 
0.770272, 0.77149, 0.771795, 0.7721, 0.772709, 0.773927, 0.774257, 
0.774588, 0.775249, 0.775579, 0.775909, 0.77657, 0.7769, 0.777231, 
0.777892, 0.779213, 0.779543, 0.779874, 0.780535, 0.780865, 0.781195, 
0.781856, 0.782186, 0.782517, 0.782847, 0.783178, 0.783508, 0.783838, 
0.784169, 0.784499, 0.78483, 0.78516, 0.78549, 0.785821, 0.786151, 
0.786481, 0.786812, 0.787142, 0.787473, 0.787803, 0.788133, 0.788464, 
0.788794, 0.789124, 0.789455, 0.789785, 0.790116, 0.790446, 0.790776, 
0.791107, 0.791437, 0.791767, 0.792428, 0.792759, 0.793089, 0.79375, 
0.79408, 0.79441, 0.795071, 0.79538, 0.795688, 0.796305, 0.796613, 
0.796921, 0.797538, 0.797846, 0.798155, 0.798771, 0.800005, 0.800313, 
0.800621, 0.801238, 0.802472, 0.804938, 0.805247, 0.805555, 0.806172, 
0.807405, 0.809872, 0.81018, 0.810488, 0.811105, 0.812339, 0.814805, 
0.815139, 0.815474, 0.816142, 0.817478, 0.820152, 0.820486, 0.82082, 
0.821488, 0.822825, 0.825498, 0.825832, 0.826166, 0.826834, 0.828171, 
0.830844, 0.83619, 0.836519, 0.836847, 0.837503, 0.838815, 0.841439, 
0.846688, 0.847016, 0.847344, 0.848, 0.849313, 0.851937, 0.857186, 
0.857492, 0.857798, 0.85841, 0.859634, 0.862082, 0.866978, 0.867284, 
0.86759, 0.868202, 0.869427, 0.871875, 0.876771, 0.877103, 0.877435, 
0.878098, 0.879426, 0.88208, 0.887389, 0.898007, 0.898317, 0.898627, 
0.899246, 0.900486, 0.902964, 0.90792, 0.917833, 0.918137, 0.918441, 
0.919048, 0.920263, 0.922692, 0.927552, 0.93727, 0.958357, 0.978034, 
0.978377, 0.978721, 0.979407, 0.98078, 0.983526, 0.989017, 0.98936, 
0.989704, 0.99039, 0.991763, 0.994509, 0.994852, 0.995195, 0.995881, 
0.997254, 0.997597, 0.997941, 0.998627, 0.99897, 0.999314, 0.999657, 
1.],
y:[6.50629e-11, 9.76969e-7, 1.95213e-6, 3.89722e-6, 7.76653e-6, 
0.0000154218, 0.0000304007, 0.0000590451, 0.000115264, 0.000160707, 
0.000161352, 0.000161995, 0.000163276, 0.000165821, 0.000170836, 
0.000180576, 0.000198904, 0.000199498, 0.000200091, 0.000201272, 
0.000203612, 0.00020821, 0.000217078, 0.000233514, 0.000233968, 
0.00023442, 0.00023532, 0.000237103, 0.000240601, 0.000247322, 
0.00024773, 0.000248136, 0.000248945, 0.000250547, 0.000253682, 
0.000259684, 0.000260078, 0.00026047, 0.00026125, 0.000262789, 
0.000265791, 0.000271489, 0.000271832, 0.000272173, 0.00027285, 
0.000274186, 0.000276783, 0.000281679, 0.000281966, 0.000282253, 
0.00028282, 0.000283938, 0.000286103, 0.000290153, 0.000290394, 
0.000290634, 0.000291108, 0.00029204, 0.000293835, 0.000297153, 
0.000297336, 0.000297517, 0.000297876, 0.000298578, 0.000299926, 
0.000300089, 0.000300251, 0.000300571, 0.000301196, 0.000302391, 
0.000302535, 0.000302677, 0.000302959, 0.000303509, 0.000304552, 
0.000304677, 0.000304801, 0.000305046, 0.000305521, 0.000306415, 
0.000306531, 0.000306645, 0.000306869, 0.000307301, 0.000308102, 
0.000308196, 0.000308289, 0.00030847, 0.000308817, 0.000309448, 
0.000309521, 0.000309592, 0.000309732, 0.000309994, 0.000310057, 
0.000310118, 0.000310237, 0.000310458, 0.00031051, 0.000310561, 
0.000310659, 0.000310706, 0.000310751, 0.000310839, 0.000310881, 
0.000310921, 0.000310999, 0.000311138, 0.000311168, 0.000311196, 
0.00031125, 0.000311275, 0.000311299, 0.000311344, 0.000311365, 
0.000311384, 0.000311403, 0.000311421, 0.000311437, 0.000311452, 
0.00031148, 0.000311492, 0.000311503, 0.000311513, 0.000311521, 
0.000311529, 0.000311536, 0.000311541, 0.000311546, 0.000311549, 
0.000311552, 0.000311553, 0.000311553, 0.000311552, 0.00031155, 
0.000311547, 0.000311543, 0.000311538, 0.000311532, 0.000311524, 
0.000311516, 0.000311506, 0.000311496, 0.000311484, 0.000311472, 
0.000311458, 0.000311443, 0.000311427, 0.00031141, 0.000311392, 
0.000311373, 0.000311332, 0.00031131, 0.000311287, 0.000311237, 
0.000311211, 0.000311184, 0.000311126, 0.000311095, 0.000311064, 
0.000310997, 0.000310852, 0.000310814, 0.000310775, 0.000310694, 
0.00031052, 0.000310474, 0.000310427, 0.00031033, 0.000310124, 
0.000310071, 0.000310016, 0.000309903, 0.000309666, 0.000309146, 
0.000309077, 0.000309007, 0.000308863, 0.000308565, 0.000307922, 
0.000307837, 0.000307752, 0.000307578, 0.000307219, 0.000306456, 
0.000306348, 0.000306239, 0.000306017, 0.000305561, 0.000304597, 
0.000302464, 0.000302322, 0.000302178, 0.000301889, 0.000301296, 
0.000300062, 0.000297398, 0.000297235, 0.00029707, 0.000296739, 
0.000296066, 0.00029468, 0.000291742, 0.000291551, 0.000291359, 
0.000290973, 0.000290192, 0.000288588, 0.000285225, 0.000284989, 
0.000284753, 0.000284278, 0.000283316, 0.000281346, 0.000277232, 
0.000268317, 0.00026803, 0.000267741, 0.000267162, 0.000265994, 
0.000263618, 0.000258714, 0.000248313, 0.000226931, 0.000201112, 
0.000174865, 0.000144445, 0.000112943, 0.0000824206, 0.0000484496, 
0.0000162373, -0.0000154971, -0.0000497517, -0.0000812325, 
-0.000114459, -0.000145811, -0.000173594, -0.00020178, -0.000225942, 
-0.000247314, -0.000247655, -0.000247995, -0.000248673, -0.00025002, 
-0.000252679, -0.000257848, -0.000267578, -0.000267849, -0.000268119, 
-0.000268658, -0.000269727, -0.000271829, -0.00027589, -0.000276138, 
-0.000276384, -0.000276876, -0.000277849, -0.000279758, -0.000283427, 
-0.000283668, -0.000283908, -0.000284386, -0.000285331, -0.000287174, 
-0.000290677, -0.000290888, -0.000291098, -0.000291514, -0.000292335, 
-0.00029393, -0.000296926, -0.000297093, -0.000297259, -0.000297589, 
-0.000298236, -0.000299488, -0.000301817, -0.000301955, -0.000302092, 
-0.000302363, -0.000302893, -0.000303909, -0.000305758, -0.000305863, 
-0.000305967, -0.000306173, -0.000306574, -0.000307331, -0.000307421, 
-0.00030751, -0.000307686, -0.000308027, -0.000308662, -0.000308737, 
-0.000308811, -0.000308957, -0.000309236, -0.000309303, -0.00030937, 
-0.000309499, -0.000309748, -0.000309807, -0.000309866, -0.00030998, 
-0.000310197, -0.000310248, -0.000310299, -0.000310398, -0.000310583, 
-0.00031063, -0.000310676, -0.000310765, -0.000310808, -0.00031085, 
-0.000310929, -0.000310967, -0.000311004, -0.000311074, -0.0003112, 
-0.000311229, -0.000311256, -0.000311307, -0.000311331, -0.000311354, 
-0.000311395, -0.000311414, -0.000311432, -0.000311448, -0.000311464, 
-0.000311478, -0.000311491, -0.000311502, -0.000311513, -0.000311522, 
-0.00031153, -0.000311537, -0.000311543, -0.000311547, -0.00031155, 
-0.000311552, -0.000311553, -0.000311553, -0.000311551, -0.000311548, 
-0.000311544, -0.000311538, -0.000311532, -0.000311524, -0.000311515, 
-0.000311505, -0.000311493, -0.00031148, -0.000311466, -0.000311451, 
-0.000311434, -0.000311398, -0.000311377, -0.000311356, -0.000311309, 
-0.000311284, -0.000311258, -0.000311201, -0.000311173, -0.000311143, 
-0.000311081, -0.000311049, -0.000311015, -0.000310944, -0.000310907, 
-0.000310869, -0.00031079, -0.000310618, -0.000310572, -0.000310525, 
-0.000310428, -0.00031022, -0.00030975, -0.000309686, -0.000309621, 
-0.000309488, -0.000309208, -0.000308593, -0.000308511, -0.000308428, 
-0.000308258, -0.000307905, -0.000307143, -0.000307034, -0.000306924, 
-0.000306699, -0.000306234, -0.000305236, -0.000305105, -0.000304972, 
-0.000304704, -0.000304149, -0.000302973, -0.00030282, -0.000302665, 
-0.000302351, -0.000301707, -0.00030035, -0.00029736, -0.000297165, 
-0.000296968, -0.00029657, -0.000295757, -0.000294063, -0.000290403, 
-0.000290162, -0.00028992, -0.000289431, -0.000288435, -0.000286374, 
-0.000281971, -0.000281703, -0.000281433, -0.00028089, -0.000279788, 
-0.000277523, -0.00027274, -0.00027243, -0.000272119, -0.000271492, 
-0.000270222, -0.00026762, -0.000262157, -0.000261774, -0.00026139, 
-0.000260616, -0.000259051, -0.000255842, -0.000249115, -0.000234402, 
-0.000233948, -0.000233491, -0.000232575, -0.000230724, -0.000226951, 
-0.000219123, -0.000202318, -0.000201779, -0.000201238, -0.000200152, 
-0.000197963, -0.000193514, -0.000184331, -0.000164815, -0.000117064, 
-0.0000656105, -0.0000646524, -0.0000636923, -0.0000617656, 
-0.000057887, -0.0000500284, -0.0000339028, -0.0000328768, 
-0.0000318486, -0.0000297858, -0.0000256345, -0.0000172284, 
-0.0000161679, -0.0000151053, -0.0000129735, -8.68382e-6, 
-7.60598e-6, -6.52596e-6, -4.35938e-6, -3.27282e-6, -2.18409e-6,
-1.09317e-6, -6.50629e-11],
mode: 'lines',
name: 'pi/8'
};

var data = [plot0,plot1,plot2];

var layout = {
  title:  'angle error',
  xaxis: { nticks: 10 },
  yaxis: { nticks: 10 },
  height: 376,
  width:  626
};

Plotly.newPlot('fig1', data, layout, {displaylogo: false, autosizable: true});

</script>


<script>
var plot0 = {
x:[2.04082e-8, 0.000306718, 0.000613415, 0.00122681, 0.0024536, 
0.00490718, 0.00981434, 0.0196287, 0.0409084, 0.0607779, 0.0802576, 
0.101388, 0.121109, 0.142481, 0.163463, 0.183035, 0.204257, 0.22407, 
0.243493, 0.264567, 0.284231, 0.305546, 0.326471, 0.345985, 0.367151, 
0.386907, 0.387242, 0.387576, 0.388245, 0.389583, 0.392259, 0.397611, 
0.408314, 0.408642, 0.408971, 0.409628, 0.410941, 0.413568, 0.418823, 
0.419151, 0.419479, 0.420136, 0.42145, 0.424077, 0.429331, 0.429638, 
0.429944, 0.430557, 0.431782, 0.434233, 0.439135, 0.439441, 0.439747, 
0.44036, 0.441586, 0.444036, 0.448938, 0.44927, 0.449602, 0.450267, 
0.451595, 0.454253, 0.454585, 0.454917, 0.455581, 0.45691, 0.459567, 
0.459899, 0.460231, 0.460896, 0.462224, 0.464882, 0.470196, 0.470506, 
0.470816, 0.471437, 0.472677, 0.475158, 0.475468, 0.475778, 0.476399, 
0.477639, 0.48012, 0.48043, 0.48074, 0.48136, 0.482601, 0.482911, 
0.483221, 0.483841, 0.485082, 0.485392, 0.485702, 0.486322, 0.487563, 
0.487873, 0.488183, 0.488803, 0.490044, 0.490348, 0.490652, 0.49126, 
0.491564, 0.491868, 0.492476, 0.49278, 0.493084, 0.493692, 0.493996, 
0.4943, 0.494908, 0.495212, 0.495516, 0.49582, 0.496125, 0.496429, 
0.496733, 0.497037, 0.497341, 0.497645, 0.497949, 0.498253, 0.498557, 
0.498861, 0.499165, 0.499469, 0.499773, 0.500077, 0.500381, 0.500685, 
0.500989, 0.501293, 0.501597, 0.501901, 0.502205, 0.502509, 0.502813, 
0.503117, 0.503421, 0.503725, 0.504029, 0.504637, 0.504941, 0.505245, 
0.505854, 0.506158, 0.506462, 0.50707, 0.507374, 0.507678, 0.508286, 
0.509502, 0.509832, 0.510162, 0.510821, 0.512141, 0.51247, 0.5128, 
0.51346, 0.514779, 0.515109, 0.515439, 0.516099, 0.517418, 0.520056, 
0.520386, 0.520716, 0.521376, 0.522695, 0.525334, 0.525664, 0.525993, 
0.526653, 0.527972, 0.530611, 0.530919, 0.531227, 0.531842, 0.533073, 
0.535536, 0.54046, 0.540768, 0.541076, 0.541692, 0.542923, 0.545385, 
0.55031, 0.550644, 0.550977, 0.551644, 0.552979, 0.555647, 0.560985, 
0.561319, 0.561652, 0.562319, 0.563654, 0.566322, 0.57166, 0.571987, 
0.572315, 0.57297, 0.57428, 0.5769, 0.58214, 0.59262, 0.592926, 
0.593231, 0.593842, 0.595064, 0.597508, 0.602395, 0.61217, 0.633371, 
0.653162, 0.672563, 0.693616, 0.713258, 0.734551, 0.754434, 0.773927, 
0.795071, 0.814805, 0.83619, 0.857186, 0.876771, 0.898007, 0.917833, 
0.93727, 0.958357, 0.978034, 0.978377, 0.978721, 0.979407, 0.98078, 
0.983526, 0.989017, 0.98936, 0.989704, 0.99039, 0.991763, 0.994509, 
0.994852, 0.995195, 0.995881, 0.997254, 0.997597, 0.997941, 0.998627, 
0.99897, 0.999314, 0.999657, 1.],
y:[-0.570796, -0.570183, -0.569569, -0.56834, -0.565877, -0.560934, 
-0.550975, -0.530769, -0.485645, -0.441914, -0.397593, -0.347977, 
-0.300337, -0.247415, -0.194362, -0.144129, -0.0891555, -0.0376999, 
0.0125049, 0.0662467, 0.115224, 0.166445, 0.214201, 0.255885, 
0.297324, 0.331864, 0.332411, 0.332957, 0.334044, 0.336204, 0.340459, 
0.348711, 0.364141, 0.364591, 0.36504, 0.365933, 0.367703, 0.371174, 
0.377839, 0.378244, 0.378646, 0.379447, 0.381031, 0.384128, 0.390034, 
0.390366, 0.390697, 0.391355, 0.392655, 0.39519, 0.4, 0.400289, 
0.400576, 0.401147, 0.402273, 0.404458, 0.40856, 0.408825, 0.409089, 
0.40961, 0.410634, 0.4126, 0.412838, 0.413075, 0.413543, 0.414459, 
0.41621, 0.416421, 0.416631, 0.417045, 0.417852, 0.419386, 0.422123, 
0.422269, 0.422413, 0.422698, 0.423249, 0.424279, 0.424401, 0.424521, 
0.424757, 0.425212, 0.426047, 0.426145, 0.426241, 0.426428, 0.426785, 
0.42687, 0.426954, 0.427117, 0.427425, 0.427498, 0.42757, 0.427708, 
0.427967, 0.428028, 0.428087, 0.428201, 0.428411, 0.428459, 0.428505, 
0.428593, 0.428635, 0.428675, 0.428751, 0.428787, 0.428821, 0.428885, 
0.428915, 0.428944, 0.428996, 0.42902, 0.429043, 0.429064, 0.429084, 
0.429102, 0.429118, 0.429133, 0.429147, 0.429159, 0.42917, 0.429179, 
0.429187, 0.429193, 0.429198, 0.429201, 0.429203, 0.429204, 0.429203, 
0.4292, 0.429196, 0.42919, 0.429183, 0.429175, 0.429165, 0.429153, 
0.42914, 0.429126, 0.42911, 0.429093, 0.429074, 0.429032, 0.429008, 
0.428984, 0.42893, 0.4289, 0.42887, 0.428804, 0.428769, 0.428732, 
0.428655, 0.428482, 0.428431, 0.428378, 0.428267, 0.428025, 0.42796, 
0.427894, 0.427755, 0.427458, 0.427379, 0.427299, 0.427133, 0.42678, 
0.425991, 0.425884, 0.425776, 0.425555, 0.425092, 0.424082, 0.423949, 
0.423813, 0.423537, 0.422964, 0.421735, 0.421585, 0.421433, 0.421125, 
0.420491, 0.419152, 0.416192, 0.415995, 0.415796, 0.415394, 0.414573, 
0.41286, 0.409158, 0.408894, 0.408628, 0.408092, 0.406999, 0.404734, 
0.399886, 0.39957, 0.399251, 0.398609, 0.397306, 0.394623, 0.388949, 
0.388588, 0.388225, 0.387495, 0.386017, 0.382988, 0.376646, 0.362852, 
0.362429, 0.362004, 0.36115, 0.359429, 0.355928, 0.348702, 0.33337, 
0.296353, 0.257633, 0.216334, 0.168413, 0.121372, 0.0684721, 
0.0178338, -0.0325052, -0.0874119, -0.138552, -0.193477, -0.24658, 
-0.295144, -0.346536, -0.393175, -0.437536, -0.484055, -0.525901, 
-0.526617, -0.527333, -0.528763, -0.531618, -0.537305, -0.548589, 
-0.549291, -0.549991, -0.551392, -0.554186, -0.559753, -0.560447, 
-0.56114, -0.562525, -0.56529, -0.56598, -0.566669, -0.568047, 
-0.568735, -0.569422, -0.57011, -0.570796],
mode: 'lines',
name: 'pi/2'
};

var plot1 = {
x:[2.04082e-8, 0.000306718, 0.000613415, 0.00122681, 0.0024536, 
0.00490718, 0.00981434, 0.0196287, 0.0409084, 0.0607779, 0.0802576, 
0.101388, 0.121109, 0.142481, 0.163463, 0.183035, 0.204257, 0.22407, 
0.243493, 0.264567, 0.284231, 0.305546, 0.326471, 0.345985, 0.367151, 
0.36746, 0.367769, 0.368386, 0.369621, 0.37209, 0.377029, 0.386907, 
0.387242, 0.387576, 0.388245, 0.389583, 0.392259, 0.397611, 0.408314, 
0.408642, 0.408971, 0.409628, 0.410941, 0.413568, 0.418823, 0.419151, 
0.419479, 0.420136, 0.42145, 0.424077, 0.429331, 0.429638, 0.429944, 
0.430557, 0.431782, 0.434233, 0.439135, 0.439441, 0.439747, 0.44036, 
0.441586, 0.444036, 0.448938, 0.44927, 0.449602, 0.450267, 0.451595, 
0.454253, 0.454585, 0.454917, 0.455581, 0.45691, 0.459567, 0.459899, 
0.460231, 0.460896, 0.462224, 0.464882, 0.470196, 0.470506, 0.470816, 
0.471437, 0.472677, 0.475158, 0.475468, 0.475778, 0.476399, 0.477639, 
0.48012, 0.48043, 0.48074, 0.48136, 0.482601, 0.482911, 0.483221, 
0.483841, 0.485082, 0.485392, 0.485702, 0.486322, 0.487563, 0.487873, 
0.488183, 0.488803, 0.490044, 0.490348, 0.490652, 0.49126, 0.491564, 
0.491868, 0.492476, 0.49278, 0.493084, 0.493692, 0.493996, 0.4943, 
0.494908, 0.495212, 0.495516, 0.49582, 0.496125, 0.496429, 0.496733, 
0.497037, 0.497341, 0.497645, 0.497949, 0.498253, 0.498557, 0.498861, 
0.499165, 0.499469, 0.499773, 0.500077, 0.500381, 0.500685, 0.500989, 
0.501293, 0.501597, 0.501901, 0.502205, 0.502509, 0.502813, 0.503117, 
0.503421, 0.503725, 0.504029, 0.504637, 0.504941, 0.505245, 0.505854, 
0.506158, 0.506462, 0.50707, 0.507374, 0.507678, 0.508286, 0.509502, 
0.509832, 0.510162, 0.510821, 0.512141, 0.51247, 0.5128, 0.51346, 
0.514779, 0.515109, 0.515439, 0.516099, 0.517418, 0.520056, 0.520386, 
0.520716, 0.521376, 0.522695, 0.525334, 0.525664, 0.525993, 0.526653, 
0.527972, 0.530611, 0.530919, 0.531227, 0.531842, 0.533073, 0.535536, 
0.54046, 0.540768, 0.541076, 0.541692, 0.542923, 0.545385, 0.55031, 
0.550644, 0.550977, 0.551644, 0.552979, 0.555647, 0.560985, 0.561319, 
0.561652, 0.562319, 0.563654, 0.566322, 0.57166, 0.571987, 0.572315, 
0.57297, 0.57428, 0.5769, 0.58214, 0.59262, 0.592926, 0.593231, 
0.593842, 0.595064, 0.597508, 0.602395, 0.61217, 0.633371, 0.653162, 
0.672563, 0.693616, 0.713258, 0.734551, 0.754434, 0.773927, 0.795071, 
0.814805, 0.83619, 0.857186, 0.876771, 0.898007, 0.917833, 0.93727, 
0.958357, 0.978034, 0.978377, 0.978721, 0.979407, 0.98078, 0.983526, 
0.989017, 0.98936, 0.989704, 0.99039, 0.991763, 0.994509, 0.994852, 
0.995195, 0.995881, 0.997254, 0.997597, 0.997941, 0.998627, 0.99897, 
0.999314, 0.999657, 1.],
y:[-0.0782914, -0.0781644, -0.0780374, -0.0777835, -0.0772761, 
-0.0762629, -0.074243, -0.0702296, -0.0616575, -0.0538284, 
-0.0463338, -0.0384252, -0.0312699, -0.0237814, -0.0167186, 
-0.0104065, -0.00388201, 0.00189119, 0.00723738, 0.0126702, 
0.0173792, 0.0220747, 0.0262554, 0.0297589, 0.033115, 0.0331605, 
0.0332058, 0.0332963, 0.0334759, 0.0338304, 0.0345199, 0.0358206, 
0.0358628, 0.0359049, 0.0359887, 0.0361549, 0.0364814, 0.0371112, 
0.038277, 0.0383108, 0.0383445, 0.0384115, 0.038544, 0.0388034, 
0.0392993, 0.0393292, 0.0393591, 0.0394185, 0.0395358, 0.0397646, 
0.0401993, 0.0402237, 0.040248, 0.0402963, 0.0403916, 0.0405771, 
0.0409281, 0.0409491, 0.0409701, 0.0410116, 0.0410935, 0.0412521, 
0.0415492, 0.0415684, 0.0415874, 0.0416251, 0.041699, 0.0418408, 
0.041858, 0.041875, 0.0419087, 0.0419747, 0.0421005, 0.0421157, 
0.0421308, 0.0421605, 0.0422184, 0.0423284, 0.0425242, 0.0425347, 
0.042545, 0.0425654, 0.0426047, 0.0426782, 0.0426869, 0.0426955, 
0.0427124, 0.0427448, 0.0428043, 0.0428113, 0.0428181, 0.0428315, 
0.0428569, 0.042863, 0.0428689, 0.0428805, 0.0429025, 0.0429077, 
0.0429128, 0.0429226, 0.042941, 0.0429454, 0.0429496, 0.0429577, 
0.0429726, 0.042976, 0.0429793, 0.0429855, 0.0429885, 0.0429914, 
0.0429968, 0.0429993, 0.0430018, 0.0430063, 0.0430085, 0.0430105, 
0.0430142, 0.0430159, 0.0430175, 0.043019, 0.0430204, 0.0430217, 
0.0430229, 0.043024, 0.0430249, 0.0430258, 0.0430266, 0.0430272, 
0.0430278, 0.0430282, 0.0430286, 0.0430288, 0.0430289, 0.043029, 
0.0430289, 0.0430287, 0.0430284, 0.043028, 0.0430275, 0.0430269, 
0.0430262, 0.0430254, 0.0430245, 0.0430234, 0.0430223, 0.0430211, 
0.0430197, 0.0430167, 0.0430151, 0.0430133, 0.0430095, 0.0430074, 
0.0430052, 0.0430005, 0.042998, 0.0429954, 0.0429899, 0.0429776, 
0.042974, 0.0429703, 0.0429624, 0.0429452, 0.0429406, 0.0429358, 
0.042926, 0.0429048, 0.0428992, 0.0428935, 0.0428816, 0.0428565, 
0.0428003, 0.0427927, 0.042785, 0.0427693, 0.0427362, 0.0426642, 
0.0426547, 0.042645, 0.0426253, 0.0425843, 0.0424966, 0.0424858, 
0.0424749, 0.0424529, 0.0424075, 0.0423116, 0.0420993, 0.0420851, 
0.0420708, 0.0420419, 0.0419828, 0.0418595, 0.0415924, 0.0415733, 
0.0415541, 0.0415154, 0.0414363, 0.0412721, 0.0409198, 0.0408968, 
0.0408736, 0.0408268, 0.0407317, 0.0405357, 0.0401197, 0.0400931, 
0.0400664, 0.0400127, 0.0399039, 0.0396804, 0.0392107, 0.0381803, 
0.0381484, 0.0381165, 0.0380523, 0.0379226, 0.0376584, 0.0371105, 
0.0359367, 0.0330377, 0.029903, 0.0264379, 0.0222506, 0.0179549, 
0.0128891, 0.00778911, 0.00245711, -0.00368126, -0.00972687, 
-0.0166044, -0.0236669, -0.0305144, -0.0382026, -0.0456096, 
-0.0530689, -0.0613649, -0.0692793, -0.0694187, -0.0695582, 
-0.0698373, -0.070396, -0.0715156, -0.0737633, -0.0739041, -0.074045, 
-0.0743269, -0.0748913, -0.076022, -0.0761635, -0.0763051, 
-0.0765883, -0.0771554, -0.0772972, -0.0774391, -0.077723, 
-0.0778651, -0.0780071, -0.0781492, -0.0782914],
mode: 'lines',
name: 'pi/4'
};

var plot2 = {
x:[2.04082e-8, 0.000306718, 0.000613415, 0.00122681, 0.0024536, 
0.00490718, 0.00981434, 0.0196287, 0.0409084, 0.0607779, 0.0802576, 
0.101388, 0.121109, 0.142481, 0.163463, 0.183035, 0.204257, 0.22407, 
0.243493, 0.264567, 0.284231, 0.305546, 0.326471, 0.345985, 0.367151, 
0.36746, 0.367769, 0.368386, 0.369621, 0.37209, 0.377029, 0.386907, 
0.387242, 0.387576, 0.388245, 0.389583, 0.392259, 0.397611, 0.397945, 
0.39828, 0.398949, 0.400287, 0.402962, 0.408314, 0.408642, 0.408971, 
0.409628, 0.410941, 0.413568, 0.418823, 0.419151, 0.419479, 0.420136, 
0.42145, 0.424077, 0.429331, 0.429638, 0.429944, 0.430557, 0.431782, 
0.434233, 0.439135, 0.439441, 0.439747, 0.44036, 0.441586, 0.444036, 
0.448938, 0.44927, 0.449602, 0.450267, 0.451595, 0.454253, 0.454585, 
0.454917, 0.455581, 0.45691, 0.459567, 0.459899, 0.460231, 0.460896, 
0.462224, 0.464882, 0.470196, 0.470506, 0.470816, 0.471437, 0.472677, 
0.475158, 0.475468, 0.475778, 0.476399, 0.477639, 0.48012, 0.48043, 
0.48074, 0.48136, 0.482601, 0.482911, 0.483221, 0.483841, 0.485082, 
0.485392, 0.485702, 0.486322, 0.487563, 0.487873, 0.488183, 0.488803, 
0.490044, 0.490348, 0.490652, 0.49126, 0.491564, 0.491868, 0.492476, 
0.49278, 0.493084, 0.493692, 0.493996, 0.4943, 0.494908, 0.495212, 
0.495516, 0.49582, 0.496125, 0.496429, 0.496733, 0.497037, 0.497341, 
0.497645, 0.497949, 0.498253, 0.498557, 0.498861, 0.499165, 0.499469, 
0.499773, 0.500077, 0.500381, 0.500685, 0.500989, 0.501293, 0.501597, 
0.501901, 0.502205, 0.502509, 0.502813, 0.503117, 0.503421, 0.503725, 
0.504029, 0.504637, 0.504941, 0.505245, 0.505854, 0.506158, 0.506462, 
0.50707, 0.507374, 0.507678, 0.508286, 0.509502, 0.509832, 0.510162, 
0.510821, 0.512141, 0.51247, 0.5128, 0.51346, 0.514779, 0.515109, 
0.515439, 0.516099, 0.517418, 0.520056, 0.520386, 0.520716, 0.521376, 
0.522695, 0.525334, 0.525664, 0.525993, 0.526653, 0.527972, 0.530611, 
0.530919, 0.531227, 0.531842, 0.533073, 0.535536, 0.535844, 0.536151, 
0.536767, 0.537998, 0.54046, 0.540768, 0.541076, 0.541692, 0.542923, 
0.545385, 0.55031, 0.550644, 0.550977, 0.551644, 0.552979, 0.555647, 
0.560985, 0.561319, 0.561652, 0.562319, 0.563654, 0.566322, 0.57166, 
0.571987, 0.572315, 0.57297, 0.57428, 0.5769, 0.58214, 0.59262, 
0.592926, 0.593231, 0.593842, 0.595064, 0.597508, 0.602395, 0.61217, 
0.612501, 0.612833, 0.613495, 0.61482, 0.617471, 0.622771, 0.633371, 
0.653162, 0.672563, 0.693616, 0.713258, 0.734551, 0.754434, 0.773927, 
0.795071, 0.814805, 0.83619, 0.857186, 0.876771, 0.898007, 0.917833, 
0.93727, 0.958357, 0.978034, 0.978377, 0.978721, 0.979407, 0.98078, 
0.983526, 0.989017, 0.98936, 0.989704, 0.99039, 0.991763, 0.994509, 
0.994852, 0.995195, 0.995881, 0.997254, 0.997597, 0.997941, 0.998627, 
0.99897, 0.999314, 0.999657, 1.],
y:[-0.0100156, -0.00999778, -0.00997993, -0.00994425, -0.009873, 
-0.00973095, -0.00944864, -0.00889123, -0.00771609, -0.00666078, 
-0.00566622, -0.00463298, -0.0037122, -0.00276252, -0.00187963, 
-0.00110091, -0.000306009, 0.000389004, 0.00102574, 0.00166612, 
0.00221583, 0.00275908, 0.00323872, 0.00363776, 0.00401756, 
0.00402269, 0.0040278, 0.004038, 0.00405826, 0.00409821, 0.00417584, 
0.004322, 0.00432674, 0.00433146, 0.00434087, 0.00435951, 0.00439613, 
0.00446668, 0.00447098, 0.00447525, 0.00448377, 0.00450062, 
0.00453367, 0.00459708, 0.00460086, 0.00460462, 0.0046121, 0.0046269, 
0.00465586, 0.00471118, 0.00471453, 0.00471786, 0.00472448, 
0.00473755, 0.00476305, 0.00481147, 0.00481418, 0.00481689, 
0.00482226, 0.00483287, 0.00485351, 0.00489254, 0.00489488, 
0.00489721, 0.00490183, 0.00491093, 0.00492856, 0.00496156, 
0.00496369, 0.0049658, 0.00496999, 0.00497819, 0.00499393, 
0.00499583, 0.00499773, 0.00500147, 0.00500878, 0.00502275, 
0.00502443, 0.0050261, 0.0050294, 0.00503583, 0.00504802, 0.00506973, 
0.00507089, 0.00507204, 0.00507429, 0.00507866, 0.0050868, 
0.00508777, 0.00508872, 0.00509058, 0.00509417, 0.00510077, 
0.00510154, 0.0051023, 0.00510378, 0.00510659, 0.00510727, 
0.00510793, 0.00510921, 0.00511164, 0.00511222, 0.00511278, 
0.00511387, 0.00511591, 0.00511639, 0.00511686, 0.00511776, 
0.00511941, 0.00511979, 0.00512015, 0.00512084, 0.00512117, 
0.00512149, 0.00512209, 0.00512237, 0.00512264, 0.00512315, 
0.00512338, 0.00512361, 0.00512402, 0.00512421, 0.00512439, 
0.00512455, 0.00512471, 0.00512485, 0.00512498, 0.0051251, 
0.00512521, 0.0051253, 0.00512539, 0.00512546, 0.00512552, 
0.00512557, 0.00512561, 0.00512564, 0.00512565, 0.00512565, 
0.00512564, 0.00512562, 0.00512559, 0.00512555, 0.00512549, 
0.00512543, 0.00512535, 0.00512526, 0.00512515, 0.00512504, 
0.00512492, 0.00512478, 0.00512463, 0.0051243, 0.00512412, 
0.00512392, 0.0051235, 0.00512327, 0.00512302, 0.00512251, 
0.00512223, 0.00512194, 0.00512133, 0.00511997, 0.00511957, 
0.00511915, 0.00511828, 0.00511637, 0.00511586, 0.00511534, 
0.00511425, 0.0051119, 0.00511128, 0.00511065, 0.00510934, 
0.00510655, 0.00510033, 0.00509949, 0.00509863, 0.00509689, 
0.00509323, 0.00508525, 0.00508419, 0.00508312, 0.00508093, 
0.00507639, 0.00506666, 0.00506547, 0.00506427, 0.00506182, 
0.00505679, 0.00504616, 0.00504478, 0.00504338, 0.00504056, 
0.00503477, 0.00502261, 0.00502104, 0.00501945, 0.00501624, 
0.00500969, 0.00499601, 0.00496636, 0.00496424, 0.0049621, 0.0049578, 
0.00494901, 0.00493078, 0.00489163, 0.00488906, 0.00488648, 
0.00488128, 0.00487071, 0.0048489, 0.0048026, 0.00479964, 0.00479667, 
0.00479069, 0.00477856, 0.00475367, 0.00470131, 0.00458627, 
0.00458271, 0.00457914, 0.00457196, 0.00455747, 0.00452793, 
0.00446661, 0.00433504, 0.00433037, 0.00432569, 0.00431628, 
0.00429731, 0.00425872, 0.00417891, 0.00400884, 0.00365412, 
0.00325957, 0.00277933, 0.0022827, 0.00169178, 0.00109108, 
0.000456717, -0.000281711, -0.00101763, -0.00186546, -0.0027481, 
-0.00361575, -0.00460414, -0.00557092, -0.0065593, -0.00767635, 
-0.00875993, -0.00877918, -0.00879844, -0.00883699, -0.00891425, 
-0.00906933, -0.00938177, -0.00940139, -0.00942103, -0.00946034, 
-0.00953911, -0.00969721, -0.00971702, -0.00973685, -0.00977654, 
-0.00985606, -0.00987596, -0.00989588, -0.00993576, -0.00995571,
-0.00997568, -0.00999566, -0.0100156],
mode: 'lines'
name: 'pi/8'
};

var data = [plot0,plot1,plot2];

var layout = {
  title:  'angular velocity error',
  xaxis: { nticks: 10 },
  yaxis: { nticks: 10 },
  height: 376,
  width:  626
};

Plotly.newPlot('fig2', data, layout, {displaylogo: false, autosizable: true});

</script>