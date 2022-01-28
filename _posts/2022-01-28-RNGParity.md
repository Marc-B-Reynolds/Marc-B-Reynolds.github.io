---
layout:       post
title:        Random number draws with fixed parity
categories:   [math]
tags:         [random,distribution]
plotly:       false
description:  Describes simple transforms to generate random numbers with a given parity and the math background.
---

\\
I tossed out a tweet about a simple method to generate a random number with either even or odd [parity](https://en.wikipedia.org/wiki/Parity_bit) (even/odd [population count](https://en.wikipedia.org/wiki/Hamming_weight)). This is a quick note for the justification.

Let's start by noting we can transform a random integer to being either even or odd by simply forcing the value of the lowest bit:


{% highlight c %}

// transform random integer 'u' to odd/even
static inline uint64_t rng_xf_odd_64(uint64_t u)  { u |=  1; return u; }
static inline uint64_t rng_xf_even_64(uint64_t u) { u &= ~1; return u; }

{% endhighlight %}

\\
These obviously return results with the desired property but there's a second consideration: we want to ensure the values are "uniform". This means that every value in the output set (the codomain) must be produced with equal probability. Well the size of both even & odd integers are half that of full set and these "maps" take exactly two input values of the domain and maps them to one value in the [codomain](https://en.wikipedia.org/wiki/Codomain) (a 2-to-1 map). So the result is uniform. A different but equally good map for the even case is:

{% highlight c %}
static inline uint64_t rng_xf_even_64(uint64_t u) { u <<= 1; return u; }
{% endhighlight %}

\\
So my tweet was basically: map the integer to odd/even and transform by the [Gray code](https://en.wikipedia.org/wiki/Gray_code) to produce odd/even parity random numbers:

{% highlight c %}

uint64_t rng_xf_parity_odd(uint64_t u)   { u  |= 1; return u ^ (u >> 1); }
uint64_t rng_xf_parity_even(uint64_t u)  { u <<= 1; return u ^ (u >> 1); }

{% endhighlight %}


\\
and the Gray code is an invertiable function (aka [bijection/one-to-one map](https://en.wikipedia.org/wiki/Bijection)) so we're good for uniformity. The final question is: why does applying the Gray code map odd/even integers to odd/even parity? There's nothing about parity on the Wikipedia page ATM. Well the "action" is really happening in the space of the inverse function. 

The inverse of the Gray code is also the binary [prefix sum](https://en.wikipedia.org/wiki/Prefix_sum). As an aside this can be computed with a carryless product. I need to complete that blog post at some point..but anyway..we don't need to actually perform a prefix sum..just look at how it behaves.

One way to think about "inverse function is the space of interest" is: We start with a uniform integer that's a prefix sum, forcing that be even/odd and then transforming to a standard binary integer. Let's take a peek at an 8-bit prefix sum:


{% highlight c %}
uint8_t prefix_sum_8(uint8_t x )
{
  x ^= x >> 4;
  x ^= x >> 2;
  x ^= x >> 1;
  return x;
}
{% endhighlight %}

\\
and if we computed `y = prefix_sum_8(x)` then the binary digits of `y` are (`+` is carryless addition, better known as `XOR`):

$$
\begin{array}{lll}
  y_7 = & x_7 \\
  y_6 = & x_6 + x_7 \\
  y_5 = & x_5 + x_6 + x_7 \\
  y_4 = & x_4 + x_5 + x_6 + x_7 \\
  y_3 = & x_3 + x_4 + x_5 + x_6 + x_7 \\
  y_2 = & x_2 + x_3 + x_4 + x_5 + x_6 + x_7 \\
  y_1 = & x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 \\
  y_0 = & x_0 + x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 \\
\end{array}
$$

\\
so each binary digit of `y` is the carryless sum of itself and all higher order digits...aka the parity.  And that's it.

Well we could wander off further into obscurity and note that we could force any bit (or set of bits) to a specific value(s) and that would constrain the given intervals (from MSB to LSB) to a given parity.  Also if we swapped all right shifts for left we'd end up with a bit-reversed gray code & suffix sum so that would flip the to the opposite direction (from LSB to MSB) and going even further we could use any invertiable carryless product and it's inverse will have some other combinations of bits...you get the rough idea...I'll stop now.


