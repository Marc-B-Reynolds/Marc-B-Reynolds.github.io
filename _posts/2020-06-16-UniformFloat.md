---
layout:       post
title:        Basic uniform random floating-point values
categories:   [math]
tags:         [random,distribution]
description:  'Micropost on unbiased uniform float generation on $\left(0,1\right]$, $\left[-1,1\right)$ and $\left[-1,1\right)$. The latter two have twice the number of samples as the standard method.'
plotly:       true
---


\\
This is *nuts-and-bolt* post. We're going to assume we have the following provided function:

* `uint64_t rng_u64()` which returns a 64-bit uniform random integer

Let's first walk though the standard method which produces uniform values on $\left[0,1\right)$. If we were to simply take our 64-bit integer and convert into a float then we'd introduce a [bias]({{site.base}}/math/2016/12/22/Pigeonhole.html). So instead we limit ourselves to the largest consecutive sequence of integers that are representable and the length of which is a power-of-two: $\left[0,2^p-1\right)$ where $p$ is the number of precision bits of the format.

{% highlight c %}
// unbiased uniform integer (in floating-point) [0,2^p-1]
inline float  rng_f32_i() { return (float) (rng_u64() >> (64-24)); }
inline double rng_f64_i() { return (double)(rng_u64() >> (64-53)); }
{% endhighlight %}

\\
and to get our desired range we simply scale the output:

{% highlight c %}
// unbiased uniform floating-point values on [0,1)
inline float  rng_f32()   { return rng_f32_i() * 0x1p-24f; }
inline double rng_f64()   { return rng_f64_i() * 0x1p-53;  }
{% endhighlight %}


\\
Which completes the *standard* method of generation (which I refer to as the *equidistant* method since that's how samples are spaced). There are a number of ways that we can flip output to $\left(0,1\right]$ but let's just show one:

{% highlight c %}
// unbiased uniform floating-point values on (0,1]
inline float  rng_f32_ez() { return fmaf(rng_f32_i(), 0x1p-24f, 0x1p-24f); }
inline double rng_f64_ez() { return fma (rng_f64_i(), 0x1p-53,  0x1p-53 ); }
{% endhighlight %}

\\
So we're simply shifting all the samples up by one position. The FMA is not required. However it keeps the operation count the same as the standard method (and the same latency if FMAs and products are the same).

When we need base uniform result on $\left[-1,1\right)$ a common method is to use the standard method to generate $u$ and then transform by $2u-1$. We can instead follow the exact same path as the methods above except start with a signed integer. This results in proceedures that have the same runtime cost and (as an added bonus) twice the number of output samples.

{% highlight c %}
// unbiased uniform signed integer (in floating-point) [-2^p,2^p-1]
inline float  rng_f32_si() { return (float) ((int64_t)rng_u64() >> (64-25)); }
inline double rng_f64_si() { return (double)((int64_t)rng_u64() >> (64-54)); }

// unbiased uniform floating-point values on [-1,1)
inline float  rng_f32_s()  { return rng_f32_si() * 0x1p-24f; }
inline double rng_f64_s()  { return rng_f64_si() * 0x1p-53;  }

// unbiased uniform floating-point values on (-1,1]
inline float  rng_f32_sez() { return fmaf(rng_f32_si(), 0x1p-24f, 0x1p-24f); }
inline double rng_f64_sez() { return fma (rng_f64_si(), 0x1p-53,  0x1p-53 ); }

{% endhighlight %}
