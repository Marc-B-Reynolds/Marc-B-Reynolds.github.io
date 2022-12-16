---
layout:       post
title:        'An evil and odious map (partition into even/odd parity in increasing order)'
categories:   [math]
tags:         [integer,parity]
plotly: false
description:  'Presents a map that paritions integers into the even and odd parity subsets where each are in increasing order'
---

\\
Integers of even and odd parity are called evil and odious respectively (first letter indicates even/odd parity...just to point out the obvious). In a previous blog post ["Random number draws with fixed parity"](https://marc-b-reynolds.github.io/math/2022/01/28/RNGParity.html) I pointed out a method of transforming even/odd (as in low bit clear/set) to even/odd parity integers. TL/DR:

* applying the [gray code](https://en.wikipedia.org/wiki/Gray_code) tranforms an even/odd integer to an even/odd parity integer
* this works in the same way as [inverse transform sampling](https://en.wikipedia.org/wiki/Inverse_transform_sampling) does (the action is in the inverse function of the one being applied)
* the inverse function of the gray code (which is also the binary prefix scan) produces the cumulative parity (from top to bottom) of input at each bit position. So the bottom bit is the parity of the input and forcing that final bit either clear or set and then transforming it by the gray code will produce an integer with the desired parity.
* I note that we can do the same thing in the opposite bit order: If the logical inverse function is the binary suffix scan then the total parity is in the top bit.
* I mentioned that the gray code and it's inverse are forms of carryless products by a constant (so if multiplying by $G$ is gray code, then the inverse gray code is multiplying by $G^{-1}$) but that all of that was in an incomplete blog post which is still *"pining for the fjords"* in my drafts directory.

\\
So we could just feed the integers $\left(n,n+1,n+2, \cdots \right)$ into the above methods and get the desired parity but it'll be missing a potentially desired property: the output won't be strictly increasing in value.  We'll visit all $b$-bit integers before moving on to all $b+1$ bit integers but within all power-of-two ranges the values will be permuted. As as example let's look at the following code snippet:

{% highlight c %}
  for(uint32_t i = 0; i</*6*/66; i += 2) {
    uint32_t evil = i ^ (i >> 1);     // gray code : cr_mul_32(i,G);
    printf("%d ", evil);
  }
  printf("\n");
{% endhighlight %}

\\
produces the following:

    0 3 6 5 12 15 10 9 24 27 30 29 20 23 18 17 48 51 54 53 60 63 58 57 40 43 46 45 36 39 34 33 96

\\
For an in-order map we can adapt a simple closed form formula by Rémy Sigrist (from April 2022) posted on OEIS (SEE: [A000069](https://oeis.org/A000069) and [A001969](https://oeis.org/A001969)):

{% highlight c %}

// odious/evil partition: partitions 32-bit integers into:
// * increasing order even parity subset OEIS: A000069 (aka evil   numbers) [0,    2^31)
// * increasing order odd  parity subset OEIS: A001969 (aka odious numbers) [2^32, 2^32)
// sign bit selects the parity, lower 31 is 'n' for the given a(n) 
uint32_t eop_code_32(uint32_t n)
{
  uint32_t r = bit_prefix_sum_32(n);   // cr_mul_32(r,-1) : inverse gray code
  return r ^ (r<<1);                   // cl_mul_32(r, 3) : bit-reflected gray code
}

// inverse function
inline uint32_t eop_code_inv_32(uint32_t n)
{
  uint32_t r = bit_suffix_sum_32(n);   // cl_mul_32(n, -1) : inverse of bit-reflected gray code
  return r ^ (r>>1);                   // cr_mul_32(n,  G) : gray code (G = bit_reverse(3))
}

// swaps the ordering to odious/evil
// * (carryless) adding one swaps the order of the two subsets
uint32_t oep_code_32(uint32_t n)
{
  return even_parity_seq_32(n) ^ 1;
}

// inverse thereof (just solve)
uint32_t oep_code_inv_32(uint32_t n)
{
  return oep_code_inv_32(n^1);
}

// OEIS: A006068 
// with carryless product ops this is: cr_mul_32(a,-1)
static inline uint32_t bit_prefix_sum_32(uint32_t a)
{
  a ^= a >>  1;
  a ^= a >>  2;
  a ^= a >>  4;
  a ^= a >>  8;
  a ^= a >> 16;
  return a;
}
{% endhighlight %}

\\
Now this is about as far as I can go without actually completing my dumb carryless product post. But anyway the sideline comments say what the logical product is at each step (these can be found in [carryless.h](https://github.com/Marc-B-Reynolds/Stand-alone-junk/blob/master/src/SFH/carryless.h)). 

The super short version is `cl_mul_32` is the standard notion of a carryless product and `cr_mul_32` is a bit reflected version.  Where 'cl' is a sum of left shifts, 'cr' is a sum of right shifts which forms a second (trivally related) commutative ring.  So we're multipling by 3 (or bit-reversed 3) and all bits set (written here as -1, but really should be `UINT32_C(~0)`) which is the multiplicative inverse of 3 for 'cl' and bit_reverse(3) in 'cr'.

I haven't thought too much about product sequences when mixed across two rings so the above might have some reductions. I also can't properly justify why it works but since the action is in the inverse function we can note why the partition occurs.  The first step of `eop_code_inv_32` is the suffix sum which results in the parity of the input being in the high bit and the second (and final) step is the gray code which leaves the high bit unmodifed.




