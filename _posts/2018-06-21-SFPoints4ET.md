---
layout:       post
title:        An empirical testing point set on sphere and caps
categories:   [math]
tags:         [testing]
description:  Code for incrementally generating a spherical Fibonacci point set on a spherical cap
plotly:       true
---

\\
I'm just tossing out some code for generating a *"Spherical Fibonacci point set"* on a spherical cap of height $h$ with some minimal notes:

* Cap is oriented around $+\mathbf{z}$
* Internal computation is incremental and performed in doubles to limit compounding of errors from skewing the distribution. The assumption is $n$ is large.
* Not generally suited for performance testing since the sequences of $z$ values is linear (assuming that matters)
* The specific SFPS variant is *NOT* skipping logical element zero which is typically the case.

<br>

{% highlight c %}
// constant turning rate: 
//   TX = cos(2pi K)
//   TY = sin(2pi K) 
//   K  = frac(phi) = 1/phi = (sqrt(5)-1)/2
const double TX = -0.73736887807831985597317725478205829858779907226562;
const double TY = -0.67549029426152362720614519275841303169727325439453;

typedef struct {
  double x,y;     // incrementally computed point on circle
  double z,dz;    // incrementally computed height on cap
} sf_walk_t;

// n = number of points to generate
// h = height of cap (ex: half-sphere=1, full-sphere=2)
void sf_walk_init(sf_walk_t* w, uint32_t n, float h)
{
  w->x  = 1.0;
  w->y  = 0.0;
  w->z  = 1.0;
  w->dz = h/n;
}

void sf_walk_next(sf_walk_t* w, float* v)
{
  double x=w->x, y=w->y;
  double ct,st;

  // current disc to cap mapping values
  ct   = w->z; 
  st   = sqrt((1.0+ct)*(1.0-ct));
  
  // output current point on cap
  v[0] = (float)(st*x);
  v[1] = (float)(st*y);
  v[2] = (float)(ct);
  
  // update point on circle: turn by 2pi*K
  w->x  = TX*x-TY*y;
  w->y  = TY*x+TX*y;

  // update height in cap position
  w->z -= w->dz;
}
{% endhighlight %}
