---
layout:       post
title:        Generating the log of a uniform random float
categories:   [math]
tags:         [random, distribution]
description:  Sketches out some possible implementations of taking the logarithm of a uniform floating point value and mentions (in passing) normal distributions.
---

\\
A common building block of sampling some random distributions is taking the [logarithm](https://en.wikipedia.org/wiki/Logarithm) of a uniform random floating-point value $u$:

$$ -\log \left(u\right)$$

\\
Let's jump back and look at how the uniform floating point value $u$ was generated:

{% highlight c %}
// double on [0,1)
inline double rng_f64(void)
{
  return (rng_u64() >> (64-53))*0x1p-53;
}
{% endhighlight %}

\\
So we generate a $b$ bit uniform integer, discard $\left(b-p\right)$ bits (where $p$ is the number of bits in the significand in the floating-point format and $b$ is the number of uniform random bits we have to work with), convert the integer to float and scale by $2^{-p}$. 

The generation method above isn't what we really want since $-\log\left(0\right) = \infty$. We could instead compute $-\log\left(1-u\right)$ but let's create a new version that yields values on $\left(0,1\right]\$:


{% highlight c %}
// double on (0,1]
inline double rng_f64_zx(void)
{
  return (1+(rng_u64() >> (64-53)))*0x1p-53;
}
{% endhighlight %}

<br>

------

Turn-key
------

\\
Now let's consider something completely different, the background of which is in a previous post [*"Higher density uniform floats"*](http://marc-b-reynolds.github.io/distribution/2017/01/17/DenseFloat.html):

{% highlight c %}

// helper functions at the end of the post

// transforms 64-bit uniform integer input 'u' into a uniform floating-point
// value on (0,1). Returns the IEEE bit pattern.
inline uint64_t pdense_bits(uint64_t u)
{
  uint32_t e = lzc64(u);  // lead zero count gives the power-of-two interval of result

  u <<= e+1;              // discard the leading zeros & the first (if any) set bit
  e   = 1022-e;           // adjust the exponent to include the binary64 exp bias
  
  // build the IEEE bit pattern
  // 1) shift the exponent into position
  // 2) place the 52 explict significand bits 
  return ((uint64_t)e<<52) | (u >> 12);
}

// transforms 64-bit uniform integer input 'u' into a uniform floating-point
// value on (0,1). Exact range is [2^-65, 1-ulp].
double pdense_01_ee(uint64_t u)
{
  return f64_from_bits(pdense_bits(u));
}

// transforms 64-bit uniform integer input 'u' into a uniform floating-point
// value on (0,1]. Exact range is [2^-65+ulp, 1].
double pdense_01_ei(uint64_t u)
{
  return f64_from_bits(pdense_bits(u)+1);
}

{% endhighlight %}

\\
We have a helper function (`pdense_bits`) which generates a partially dense uniform floating point value on $\left(0,1\right)$ (as a IEEE bit pattern) and a pair of user functions.  The first simply converts the output to a float and the second shifts by 1 ULP to give the output range $\left(0,1\right]$. Let's spitball where we're at:

* standard method: generates all representiable FP values on $\left[1/2,1\right]$ and each successively smaller POT interval has half as many samples as the previous. The total number of samples generated is $2^{53}$
* partially dense method: generates all representiable FP values on $\left[2^{-12},1\right]$ and each successively smaller POT interval has half as many samples as the previous. The total number of samples generated is $2^{12}\left(2^{53}+1\right)+1$

\\
If we were to transform these by $-\log_2$ then the output ranges are $\left[0,53\right]$ and $\left[0,65\right]$ respectively. So we explode the total number of samples and extend the tail by about 20.75%.  To eyeball the cost of the two transforms see here on [godbolt](https://gcc.godbolt.org/z/5VtdA3).

<br>

------

Some assembly required <small> (really quite alot)</small>
------

\\
In this section I'm going to sketch some things we *could* do to tweak performance and/or accuracy of the previous. I'm not going to carry though any of it..'cause computing logs in software is a PITA and I'm not really motivated to talk about NA and function approximations.

Let's consider some basic [log identities](https://en.wikipedia.org/wiki/List_of_logarithmic_identities). First changing the base to base-2:


$$
\log \left(x\right) = \log \left(2\right) \log_2 \left(x\right)
$$

\\
and from the product:

$$
\log_2 \left(2^e~i\right) = e+ \log_2 \left(i\right)
$$

\\
Now let's rip apart our helper function from the previous section:

{% highlight c %}

double rng_log_xform(uint64_t u)
{
  // exponent for values on [1/2, 1)
  const uint64_t half_bits = UINT64_C(0x3fe0000000000000);
  
  uint32_t e = lzc64(u);
  u <<= e+1;                      // eat-up leading zeroes and first set bit (if any)
  u = half_bits | (u >> 12);      // put significand in place and OR in the exponent

  // d on (1/2, 1]
  double d = f64_from_bits(u+1);  // lose the +1 for [1/2,1)

  // we can use either natural or base-2 log
#if defined(USE_LOG2)
  return log_2*(e-log2(d));       // could use Kahan's ab+dc method
#else
  return log_2*e-log(d);          // fma(log_2,e,-log(d));
#endif
}

{% endhighlight %}

\\
This is just a sketch and I'm only gonna do some basic bullet points. 

* We have partially reduced input for log (and we know it's not a denormal or a special) so we can chop off some low hanging fruit.
* We're integer up to needing to compute the log. This kinda interesting if we've got an accuracy jonze going on. A very popular thing to do is to use Tang's method. This starts by examining the upper bits:  [*musl*](https://musl.libc.org/) does this (SEE: [log2.c](https://git.musl-libc.org/cgit/musl/tree/src/math/log2.c) & [log2_data.c](https://git.musl-libc.org/cgit/musl/tree/src/math/log2_data.c)) and it's talked about in [*"Computing floating-point logarithms with fixed-point operations"*](https://hal.inria.fr/hal-01227877/).

There's a just ton of moving parts here unless you're just going to spitball the result. What's coming next should be a really big one. Local examples include: which base for the log? What error bound/cost target? Single value or unevaluated-pair result from it? How do we combine all these partial results together?

<br>

------

Helper functions
------


{% highlight c %}

inline double f64_from_bits(uint64_t x)
{
  double d; memcpy(&d,&x,8); return d;
}


// lead zero count wrapper (this isn't my fault).
// 1) lead zero counting that handles (in a defined manner) a zero input
//    was introduced in Haswell.
// 2) __builtin_clzl was define in terms of the pre Haswell instruction
//    so is treated as undefined for zero input (even if compiling for
//    Haswell+)
// 3) clang does some aggressive constant folding and if see a const zero
//    input it will generate whacked output. BUT it also (seems) to always
//    observe and remove the unneeded "test" below on haswell+.
// 4) GCC doesn't remove the unneeded "test" and doesn't seem to generate
//    hacked constants like clang. This will surely break someday so the
//    safer thing to do is change this to same as clang and hope a compiler
//    fix comes soonish.
// 5) MS just uses intel intrinsics. Directing all post haswell intel
//    to this is probably the best option.
#if defined(__clang__)
inline int lzc64(uint64_t x) { return x ? __builtin_clzl(x) : 64; }
#elif defined(__GNUC__)
inline int lzc64(uint64_t x) { return __builtin_clzl(x); }
#else
#include <immintrin.h>
inline int lzc64(uint64_t x) { return _lzcnt_u64(x); }
#endif
{% endhighlight %}




