---
layout:       post
title:        Uniform points in sphere and capped cone
tagline:      
categories:   [distribution]
tags:         [distribution, random]
description:  micro post on title
plotly:       true
---

In a recent blog post[^post] (inspired from reading Peter Shirley's *Ray Tracing in One Weekend*) Karthik Karanth derives generating uniform points inside the unit sphere from scratch.  This is just a quick note on how this can be improved.


------

Uniform points in sphere
------

\\
If we're given a black-box implementation of point on unit sphere we could adapt that to an interior point by properly scaling the result.  We need a uniform parameterization of the volume: the volume of a sphere is $\frac{4}{3}\pi~r^3$, so the relative cumulative volume is $r^3$ which give a scale factor which is the cube-root of a uniform float.

In a previous post[^uniform] I sketched out generating uniform points on sphere as transformed points in disc.  Adding the scaling factor to the example code gives:

{% highlight c %}
void uniform_s2_filled(vec3_t* p)
{
  float d,s,r;
  vec2_t v;

  d    = uniform_disc(&v);      // point on disc
  r    = cbrtf(rng_f32());      // uniform param of radius
  s    = 2.f*sqrtf(1.f-d)*r;    // composed scales for x,y
  p->x = s*v.x;
  p->y = s*v.y;
  p->z = r*(1.f-2.f*d);
}
{% endhighlight %}

\\
The `uniform_disc` function from the previous post is a rejection method with rejection rate of ~0.21.  A drop in replacement without rejection would look like:

{% highlight c %}
float uniform_disc(vec2_t* p)
{
  float d = rng_f32();
  float a = 2.f*PI*rng_f32();
  float r = sqrtf(d);
  p->x = r*cosf(a);
  p->y = r*sinf(a);

  return d;                   // r^2
}
{% endhighlight %}

\\
So a disc transform based method drops the inverse trig op, one pair for forward trig ops without rejection and is trig free for rejection based.

<br>

------

Uniform points in spherically capped cone
------

\\
Tossing out another adaption of my previous post:  If we take the uniform point on spherical cap and rescale as above we get a uniform point in the volume of the intersection of the unit sphere and conic:


{% highlight c %}
// 'h' is height of the cap
void uniform_intersect_s2_cone(vec3_t* p, float h)
{
  float k,s,r;
  vec2_t v;

  r = cbrtf(rng_f32());
  k = h*uniform_disc(&v);
  s = sqrtf(h*(2.f-k))*r;
  p->x = s*v.x;
  p->y = s*v.y;
  p->z = r*(1.f-k);
}
{% endhighlight %}


------

References and Footnotes
------

[^post]:    *Generating Random Points in a Sphere*, Karthik Karanth, 2018 [link](http://karthikkaranth.me/blog/generating-random-points-in-a-sphere/)
[^uniform]: *Uniform points on disc, circle, sphere and caps*, 2016 [local post](http://marc-b-reynolds.github.io/distribution/2016/11/28/Uniform.html)

