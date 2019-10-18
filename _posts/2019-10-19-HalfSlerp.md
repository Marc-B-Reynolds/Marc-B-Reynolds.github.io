---
layout:       post
title:        "A slerp by half"
categories:   [quaternions]
tags:         [interpolation]
description:  'Some notes computing the rotation half away between two inputs'
plotly: false
---

\\
Okay this is a very sloppy thing...it all started when this tweet came across my timeline:

<blockquote class="twitter-tweet"><p lang="en" dir="ltr">Did you know all rigid transformations in 2D can be parameterized as a rotation about a pivot? Solving for the pivot point of a transformation is a fun brain-teaser you can inflict on your friends!<a href="https://t.co/XkppNf2cuV">https://t.co/XkppNf2cuV</a> <a href="https://t.co/nsFDj3YweZ">pic.twitter.com/nsFDj3YweZ</a></p>&mdash; Johnathon Selstad (@JohnSelstad) <a href="https://twitter.com/JohnSelstad/status/1184136857598087168?ref_src=twsrc%5Etfw">October 15, 2019</a></blockquote> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>

\\
which I found pretty neat so I took a quick peek at the source code. Now, me being me, I tossed out an observation which will not be generally obvious.  (Note it doesn't matter in this case...it's not a computation that's going to be statistically significant nor is tighten the error bound going to be of any interest...the point of the exercise it so demo a property.) Specifically one operation performed is a *SLERP* with a fixed interpolation value of $\frac{1}{2}$.  Basically I made the "clever" observation that algebraically we have:

$$ \text{slerp}\left(A,B,t\right) = \left(BA^*\right)^t~A  $$

\\
and since $t=\frac{1}{2}$ this reduces to:

$$ \left(BA^*\right)^{1/2}~A  = \sqrt{BA^*}A $$

\\
The code also fixes $A=1$ so this further reduces to $\sqrt{B}$. Using trig half-angle identities we get:

{% highlight c %}
// square root of unit quaterion Q. If halving rotation angle also requires w >= 0
//   Q = cos(t)+sin(t)U, |Q|=1
//   sqrt(Q) = cos(t/2) + sin(t/2)U
static inline void quat_usqrt(quat_t* q)
{
  float d = 1.f + q->w;       // 1+cos(t)
  float s = 1.f/sqrtf(d+d);   // 1/sqrt(2+2cos(t))
  q->w = d*s;                 // cos(t/2) = (1+cos(t))/sqrt(2+2cos(t))
  quat_bv_scale(q, s);        // {x,y,z} *= s, s = sin(t/2)/sin(t)
}
{% endhighlight %}

\\
As noted in the code comments this version does not handle input with negative $w$. This can be patched up by taking the abs of $w$ when computing $d$ and scaling the bivector by the sign of $d$ (see below).

I quicky remembered that most people won't be familiar with quaternion square roots, so doubled-down with a follow-up "clever" comment along the lines of: *"Oh and if we think of our two quaternions as points on the 3-sphere then they along with the origin form a plane...so just bisect the coord!"*.  Now of course by "clever" I mean the formal meaning of completely useless to anyone expect those that already know what's being said.


------

Being not clever <small>and unfix A</small>
------

\\
What I really should have said is that *nlerp* is exact at $t=\frac{1}{2}$ so one could just use that instead of slerp. That works out to be about the same complexity as the above (assuming *nlerp* is inlined) with the added patch of handling input with negative $w$.  But since I'm typing let's sketch out the second comment...with a busy figure:


![rot]({{site.base}}/assets/figures/misc/qbisect.svg 'go figure!'){: .center-image }.

Super quick:
* We have two rotations represented by unit quaternions $A$ (blue points) and $B$ (red points).
* There are two of each since all points on a line through the origin (excluded) represent the same rotation. Since we're restricting ourselves to unit magnitude we're down to two choices for a given rotation.
* $A$, $B$ and zero (origin) form a plane. This plane intersected with the set of unit quaternion forms a unit circle. We're a 2D problem from here. 
* LERPing from $A$ to $B$ linearly parameterizes the purple coord that connects them. The small purple point is at 1/2.
* SLERPing from $A$ to $B$ linearly parameterizes the angle of the arc between them. The large purple point is at 1/2.
* The line through the origin and $A+B$ (other small purple point) bisects both the coord and the arc.

We have two remaining gotchas:
* We need to choose the shortest path between the pair. If $ A\cdot B $ is negative we need to negate one of the inputs. Gray coord show the long path.
* When the angle between them is $\pi/2$, when \left(A \cdot B  = 0 \right)$, we're at a branch point, the two possible paths are the same length and we have to choose one. This isn't a quaternion specific issue.

\\
So we need to compute $\left(A+B\right)~$ if $A \cdot B \ge 0$ and otherwise $\left(A-B\right)$, then normalize that result:

{% highlight c %}
// compute normalize(A+B) taking "shortest path" into account
void quat_bisect(quat_t* r, quat_t* a, quat_t* b)
{
  float d = quat_dot(a,b);
  quat_wsum(r, a, b, 1.f, sgn(d));  // r = A + sign(A.B)B
  quat_normalize(r);
}
{% endhighlight %}


\\
We can rework the previous by noting the normalization scale is: $1/\sqrt{d}$ with

$$ d = (A+B)\cdot(A+B) = A \cdot A + B \cdot B + 2~A \cdot B = 2 + 2 ~A \cdot B $$

{% highlight c %}

void quat_bisect(quat_t* r, quat_t* a, quat_t* b)
{
  float d  = quat_dot(a,b);
  float sa = 1.0f/sqrtf(2.f+2.f*fabsf(d));
  float sb = sgn(d)*sa;
  quat_wsum(r,a,b,sa,sb);
}

{% endhighlight %}

\\
These two variants have approximately the same number and mix of operations but with the shuffling around dependencies. Just for fun we could pull out the scaling by the normalization factor and return it to be folded in a later operation:


{% highlight c %}

// R = A+B, returns 2(1+A.B), taking shortest path into account
float quat_bisect_x(quat_t* r, quat_t* a, quat_t* b)
{
  float d = quat_dot(a,b);
  quat_wsum(r, a, b, 1.f, sgn(d));  // r = A + sign(A.B)B
  return 2.f + fabsf(d+d);
}

{% endhighlight %}


------

The helper functions
------

<br>

{% highlight c %}

// return x >= 0 ? 1.f : -1.f
inline float sgn(float x) { return copysignf(1.f,x); }

typedef union { struct{ float x,y,z,w; }; float f[4]; } quat_t;

inline void quat_set(quat_t* q, float x, float y, float z, float w)
{
  q->x=x; q->y=y; q->z=z; q->w=w;
}

// scale bivector part by 's': Q = (a+B) -> (a+sB)
inline void quat_bv_scale(quat_t* q, float s)
{
  q->x *= s; q->y *= s; q->z *= s;
}

// sQ
inline void quat_scale(quat_t* q, float s)
{
  q->x *= s; q->y *= s; q->z *= s; q->w *= s;
}

// A.B
inline float quat_dot(quat_t* a, quat_t* b)
{
  return a->x*b->x + a->y*b->y + a->z*b->z + a->w*b->w;
}

// A.A
inline float quat_norm(quat_t* a) { return quat_dot(a,a); }

// weighted sum:  R = sa A + ab B
inline void quat_wsum(quat_t* r, quat_t* a, quat_t* b, float sa, float sb)
{
  float x = sa*a->x + sb*b->x;
  float y = sa*a->y + sb*b->y;
  float z = sa*a->z + sb*b->z;
  float w = sa*a->w + sb*b->w;
  quat_set(r,x,y,z,w);
}

inline void quat_normalize(quat_t* r)
{
  float s = 1.f/sqrtf(quat_norm(r));
  quat_scale(r,s);
}

{% endhighlight %}

<br>


