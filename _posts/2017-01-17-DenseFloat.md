---
layout:       post
title:        Higher density uniform floats
tagline:      
categories:   [distribution]
tags:         [random]
description:  a brief note on generating uniform floats
---

------


Prelim
------

We will start with the following assumptions.  We have a pseudo random number generator with a period of $2^P$ which returns a uniform (in all bits) bit-width result per state update.

For a 32-bit generator, the basic result returns a uniform integer: `uint32_t rng_u32()`

Other bit-widths follow naturally assuming that they are the full-width of the working register. Also I will only talk about singles.  A simple replacement of appropriate values/terms extends to doubles (or other power-of-two FP representations).  I will assume familarity with IEEE binary formats.

<br>


------

Equidistant uniform <small>what all generators do</small>
------

\\
The standard method of generating uniform floats:

{% highlight c %}
float rng_f32(uint32_t n) { return (rng_u32() >> 8) * 0x1p-24f; }
{% endhighlight %}

\\
Breaking this down we start with a 32-bit integer.  Shifting right[^mask] by eight leaves us with 24 unbiased and uniform bits and 24 is the number of bits of precision a single can hold in the normal range. Finally we multiply by $2^{-24}$ giving us with a uniform (and unbiased) single on $\left[0,1\right)$

Another way to look at it is we are producing $2^{24}$ distinct values of the form $\frac{u}{2^{24}}$ with integer $u \in \left[0, 2^{24}-1\right]$.  These values are all *equidistantly* spaced across the range.  This is the standard method since it's both cheap and commonly what one needs.

<br>

#### Break-it down

If we think about how binary floating point values are represented, we have a sign bit (we'll ignore), a power-of-two exponent, and a significand (including the implied bit).

* A result in top half of the range $\left[.5, 1\right)$ occurs with a probability of .5, the exponent (without offset[^offset]) is -1 and the *equidistant* method geneates all $2^{23}$ representable values on the range.  We are in this range if the top bit of $u$ is set.

* If the top bit of $u$ is clear but the second bit is set, then were on $\left[.25, .5\right)$, the exponent is -2 and it generates half the representable values on the range.

* If the top two bits are clear but the third is set, then were on $\left[.125, .25\right)$, the exponent is -3 and it generates one quarter of the representable values on the range.

* At the other end we have one value with exponent of -24, two with -23, four with -22, etc. (setting aside the $u=0$ case)

This all boils down to:

* We have uniform bits
* The probability any given bit is clear is .5, so the probability that the top $n$ bits[^anybits] are clear is $2^{-n}$.
* This is the geometric distribution[^gdist] with $p=.5$
* The conversion of integer to float implictly performs a lead zero count
* We're dropping one bit of density per exponent because of the lead zero count/geometric distribution.

Of course none of this matters unless you actually need uniform results with higher density.

<br>

------

Dense uniform floats
------

\\
Defining an toy problem statement: Generate a uniform float on $\left[0,1\right)$ which produces all representable values on the interval. 

<div class="alert alert-danger" role="alert" markdown="1">
All the code snippets are example toys. The final section is important.
</div>

Instead of throwing away bits and then using the automatic lead-zero-count of integer to float conversion we can manually perform the LZC to generate dense uniform floats.  The *new-jack* CS way is to break the problem down into easy/reusable parts.  So first we need to generate a random value with a geometric distribution with .5 probability of success:

{% highlight c %}
uint32_t rng_geometric()
{
  // structured to work with intel bsr w/o extra ops
  uint32_t c = 0;
  uint32_t u;
  do { u = rng_u32(); c++; } while(u==0);
  c--;
  return 33*c + lzc(u);
}
{% endhighlight %}

\\
The loop is ensured to terminate. On any given iteration the probabilty we fall through is $1-2^{-32} \approx 0.9999999998 $. Note that LZC on a 32-bit integer returns a value on $[0,32]$ which is 33 values..an easy thing to get tripped up on.

We can rewrite[^bprob] to eliminate uneeded overhead for the statistically significant case:

{% highlight c %}
uint32_t rng_geometric()
{
  // still okay for bsr w/o extra ops
  uint32_t u = rng_u32();
  
  if (u != 0) { return lzc(u); }

  // we reach here ~2.3e-8 percent of the time. If we
  // switched to a 64-bit generator it'd be ~5.4e-18%
  uint32_t c = 0;
  do { u = rng_u32(); c++; } while(u==0);
  return 33*c+lzc(u);
}
{% endhighlight %}

\\
Taking the definition of binary32 format our dense float routine becomes:

{% highlight c %}
// uniform on [0,1)
float rng_toy_1()
{
  uint32_t e = rng_geometric();            // compute the exponent

  // normal single exp on [-1,-126]
  if (e <= 125) {
    uint32_t s = rng_u32() >> 9;           // 23 bit for significand
    return f32_from_bits((126-e)<<23|s);   // convert bit representation to fp
  }
 
  // the probabilty we reach here is 2^-126 or 1.2e-36%

  // denormal (exp with offset is zero)
  if (e <= 148) {
    uint32_t s = rng_u32() >> (e-(126-9)); // no implied bit, shift way extra exp
    return f32_from_bits(s);
  }

  // we get zero 1.4e-43% of the time
  return 0.f;
}
{% endhighlight %}

\\
Which isn't too bad, not much over embellishment. Just needs some comments and meanful numeric defines. 

Let's do something similar in an *OG* style.  Rather than recreate the same thing, let's generate a value on $\left(0,1\right)$ with all representable value on the normal range (the smallest value will be $2^{-126}$).  Not only are denormals highly unlikely to occur they are probably useless or undesirable results.

First we'll generate one uniform integer and try to get as far as we can with it. We need 23 bits for the significand, so that leave us 9 to try to compute the exponent. If that didn't work out can we complete with exactly one more?  And if THAT doesn't work out, well we don't care very much.

{% highlight c %}
// uniform on (0,1)
float rng_toy_2()
{
  // not bsr friendly as written..choosing clarity
  uint32_t u = lzc(rng_u32());
  uint32_t e,s;
  
  if (e <= 8) {
    e  = 126-u;
    u &= 0x7FFFFF;  // 23 bit for significand
    return f32_from_bits(e<<23|u);
  }

  // here:   2^-9  or ~0.195313%
  // 64-bit: 2^-40 or ~9.1e-11% 

  // if we have the exponent we can complete
  if (u != 0) {
    e = 126-u;
    s = rng_u32() >> 9;
    return f32_from_bits((126-e)<<23|s);
  }

  // here:   2^-32 or ~2.33e-8%
  // 64-bit: 2^-64 or ~5.42e-18% 
  e = 33 + rng_geometric();
  if (e > 148) e = 148;       // cut the distribution

  s = rng_u32() >> (e-(126-9));
  return f32_from_bits(s);
}
{% endhighlight %}

\\
Statistically we done after generating one uniform instead of the minimum two in the previous.  There are many possible variants possible here.

<br>

------

The parts I'm not tell you
------

\\
Aside: Although written in terms of leading zeroes any choice of leading/trailing ones/zeroes is okay.

As written the results will be biased[^bias].  Not due to the pigeonhole principle, but the because of the period of the generator.  A given assumption is that we have a generator with period $2^P$.  Considering the *state* of the generator as a $P$ bit integer, then the sequence produced is a permutation of $\left[0,2^P-1\right]$. Specifically each value will occur exactly once in the sequence.  As an example let's say we're using a 32-bit generator with $P=32$ (a single state chunk).  Zero be returned exactly once per period by `rng_u32()`, so our geometric generator would return exactly one value greater than 31 (based on whatever number follows zero in the sequence) once per period and never any greater values.  Additionally many current generators have periods of $2^P-1$, where zero is excluded from the sequence.  In the $2^{32}-1$ case the geometric generator would never return a value greater than 31. For an output value of $n$ there would be $2^{31-n}$ positions in the sequence that generate that value.  The few the positions the less random these results will be.

Possible solutions:

* Don't care and structure the code to not deal with impossible situations.
* Use a generator with larger $P$.
* Don't use all the returned bits, such as XOR the top half with bottom half and use that half bit-width result.
* Extension of previous: generate two values, XOR and use that result instead.


<br>

------

References and Footnotes
------

[^offset]:  Using *offset* instead of *bias* to prevent confusion.
[^bias]:    I don't know how to measure this bias. :(
[^bprob]:   If language/compiler supports probability hints then adding the hint is reasonable.
[^mask]:    Under the assumption of uniform in all bits, masking off the low bits is equivalent.
[^anybits]: Likewise any for any choice of $n$ bits being all cleared or set.
[^gdist]:   *"Geometric distribution"*, Wikipedia. ([link](https://en.wikipedia.org/wiki/Geometric_distribution))
