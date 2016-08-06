---
layout:       post
title:        Uniform points on disc, circle, sphere and spherical cap
tagline:      
categories:   [distribution, random]
tags:         [foundation]
description:  a brief on some uniform spatial distributions
---

Here I will only consider pseudo-random point generation.  Well-spaced point generation is a different topic.  The code examples are intended for clarity.

------

Uniform points on the unit disk
------

\\
One method to generate a uniform random point $p$ on the unit disk ($\mathbb{D}^1$) is as follows.  Generate two uniform values: $r \in \left[ 0,1 \right] $ , $ \theta \in \left[ -\pi, \pi \right) $ and apply[^mwDisc]:

$$ p = \left(\sqrt{r}\cos\left(\theta\right), ~\sqrt{r}\sin\left(\theta\right)\right) $$

\\
Using symmetry and trig identities we can drop this down to one trig approximation and pair of square roots.  Instead I will only talk about using the rejection method.

<br>

#### Rejection method

Choose a 2D point $p$ uniformly in each dimension: $p_x,p_y \in \left[ -1,1 \right] $ and reject values when outside the unit disc.  The average number of iterations is $\frac{4}{\pi} \approx 1.27324$.  An example implemenation that returns $p \cdot p$ as well as the generated point:

{% highlight c %}
float uniform_disk(vec2_t* p)
{
  float d,x,y;

  do {
    x = 2.f*next_f32()-1.f;   // (-1,1)
    y = 2.f*next_f32()-1.f;   // (-1,1)
    d = x*x + y*y;            // p.p
  } while(d >= 1.f);

  p->x = x;
  p->y = y;

  return d;
}
{% endhighlight %}


<br>

------

Uniform points on unit circle
------

\\
Points on the unit circle 

John von Neumann presented a techinque to transform uniform points on the unit disc to uniform points on the circle in 1951[^jvnCircle].  If we have a uniform point $p = \left(x,y \right)$ on the unit disc, then:

$$ \frac{p^2}{p \cdot p} = \frac{1}{p \cdot p}  \left(x^2-y^2, 2xy \right)$$


<br>

------

Uniform points on the unit sphere 
------

\\
In 1972 George Marsaglia[^gm] introduced the $\left(\mathbb{S}^2\right)$

Generate a uniform point on the unit disk $\left(x,~y\right)$


$$
   \left(2x\sqrt{1-\left(x^2+y^2\right)},~2y\sqrt{1-\left(x^2+y^2\right)},~1-2\left(x^2+y^2\right)\right)
$$


{% highlight c %}
void uniform_s2(vec3_t* p)
{
  float d,s;
  vec2_t v;

  d = uniform_disk(&v);  
  s = 2.f*sqrtf(1.f-d);
  p->x = s*v.x;
  p->y = s*v.y;
  p->z = 1.f-2.f*d;
}
{% endhighlight %}



<br>

------

Uniform points on the annulus
------

<br>

------

References and Footnotes
------

[^mwDisc]:   *"Disk Point Picking"*, Eric W. Weisstein, [MathWorld](http://mathworld.wolfram.com/DiskPointPicking.html)
[^mwCircle]: *"Circle Point Picking"*, Eric W. Weisstein, [MathWorld](http://mathworld.wolfram.com/CirclePointPicking.html)
[^jvnCircle]: *"Various Techniques Used in Connection With Random Digits"*, John von Neumann, 1951 ([PDF](https://dornsifecms.usc.edu/assets/sites/520/docs/VonNeumann-ams12p36-38.pdf))
[^gm]:     *"Choosing a Point from the Surface of a Sphere"*, George Marsaglia, 1972 ([PDF](http://projecteuclid.org/euclid.aoms/1177692644))