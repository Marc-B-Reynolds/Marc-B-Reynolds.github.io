---
layout:       post
title:        Quick-n-dirty permutations
tagline:      for empirical testing and profit
categories:   [math]
tags:         [random,sequence,permutation]
description:  Some easy methods for choosing unique integers
---

Being able to draw $n$ unique integers from a set of $m$ is a handy tool-box function to have available.  Some number of times I've commented on a blog post or twitter: "You could also just use an RNG here".  Without some supporting details comments like this frankly are not very helpful so this post intended to cover the brushstrokes.

The subject is not "perfect shuffle" such as produced by *Fisher-Yates shuffle*[^fisheryates] variants.  The intended cases are when it would be overkill or undesireable to impossible to perform (say due to the size of $n$).

------


Problem statement <small>very roughly</small>
------

\\
Given the set $M=\left[0,~m\right)$ we want to form a periodic sequence(s) $S$ with period $m$ such that each element of $M$ appears exactly once...oh and it'd be rather nice if the sequence is rather random-ish.

<br>

------

Simple case <small>use an off the shelf PRNG</small>
------

\\
The simplest and probably most interesting case is when $m$ is (approximately) register size so say any 32 or 64 bit integer value.  I'll stick to use 32-bit for examples.  In that case we can simply use a fair number of non-cryptographic pseudo random number generators.  Let's just jump in with some strawman code:

{% highlight c %}

typedef struct { uint32_t state; } lsc_state_t;

static inline uint32_t next(lsc_state_t* s) 
{ 
  uint32_t state = s->state;       // previous state
  s->state = LCG_A*state + LCG_B;  // full-period bijective state update
  return state;                    // return previous state
}

{% endhighlight %}


\\
This is a classic 32-bit linear congruential generator[^lcg] (or sequence), although I've slightly restructured how `next` behaves.  Specifically it's returning the value of `state` prior it being updated.  Setting that aside we have 32-bits of state, which is updated by a special form of a bijective function (one that returns a full-period sequence) and it returns the full state.  This gives a scheme to grab $n$ unique integers from $\left[0,2^{32}\right)$.  There are a fair number of base methods we could use instead but I'm going to mostly focus on LCGs. LCGs have a quite bad reputation so I'll toss out a few point.  Part of their bad rep is due to historic bad choices of multiplicative constants which isn't an issue today since we have freely available testing tools or better yet can simply grab the values from a table[^tables]. Peter Grogono wrote about his poor constant choice and how it came to be used in a blog post[^grogono] and there's always RANDU.  Another problem is "Falls mainly in the Planes" where if we take $d$ sequential elements as treat them as points in d-dimensional cube[^lcs] then all the points fall on a finite number of planes. Not an issue here since our problem is one dimensional.




\\
xxx

### Configurable permuations

xxx

<br>

------

Example 1
------

<br>

------


References and Footnotes
------

[^lcs]:  *"Linear congruential sequences"*, Local post ([page]({{site.base}}/math/2017/09/12/LCS.html))
[^lcg]:  *"Linear congruential generators"*, Wikipedia ([page](https://en.wikipedia.org/wiki/Linear_congruential_generator))

[^fisheryates]: *"Fisher-Yates shuffle"*, Wikipedia ([page](http://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle))

[^grogono]:      *"The Grogono Generator: an Apologia"*, Peter Grogono 2013 ([page](http://users.encs.concordia.ca/~grogono/RNG/grog-gen.html))

[^tables]:       *"Tables of linear congruential generators of different sizes and good lattice structure"*, Pierre L'Ecuyer, 1999 ([PDF](http://www.ams.org/journals/mcom/1999-68-225/S0025-5718-99-00996-5/S0025-5718-99-00996-5.pdf))

[^randu]: *"RANDU"*, Wikipedia ([page](https://en.wikipedia.org/wiki/RANDU))











