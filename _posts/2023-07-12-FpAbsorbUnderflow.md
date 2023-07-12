---
layout:       post
title:        'A fast path floating point underflow test elision from absorption'
categories:   [math]
tags:         [float]
description:  'A sketch of a fast path filter to avoid explicit underflow checking following an addition or subtraction.'
plotly: true
---

\\
Okay. Sorry for the click-bait title. The contents aren't going be as exciting as you might expect. We're going to sketch out a strawman problem where:

1. We need to detect and handle underflow to meet accuracy requirements for the given computation
2. We need a high performance implementation and are willing to potentially sacrifice KISS (YMMV)

It's going to be *very* strawman because the details of the actual computation and what we know are rather quite important. The problem statement is will be stripped down to the bare minimum to demonstrate the base idea.  The sketch in code:

{% highlight c %}

  // optimistically perform a series of computations
  // that compute 'a' & 'b' assuming no input specials,
  // no overflow occurs and no underflow occurs that
  // can change the result of the follow addition or
  // subtraction.
  //
  // REQUIREMENT: if either a or b has potentially lost
  // precision in its computation then that value MUST
  // at this point be either a denormal or zero.
  
  // one or the other of these (called that r)
  s = a + b;
  d = a - b;

  // assuming r can be positive or negative and both
  // signs need to fall into the same common case
  t = fabs(r);

  // assuming we're near the top of the function we
  // can capture special inputs and hitting overflow
  // in a single high probability test. (assuming all
  // inputs are used in a & b and stuff like no divide
  // is mapping an infinity to zero..and it matters)
  if (t < infinity) {

    // if 't' is greater than the "magic" number
    // threshold then it's impossible that it
    // has lost any precision due to one of the
    // inputs being a denormal.  The magic number
    // is very quite small.
    if (t > threshold) {
    }

    // we can't determine if 'r' or |r| could have 
    // lost precision by only inspecting it. We
    // need to dig deeper to determine.
  }

  // any handling of specials and/or overflow put here
{% endhighlight %}

\\
So the remainder of this post is come up with the *magic* number `threshold` and giving a rough ideal how efficient this kind of fast-path testing strategy might be.



<br>

------

From floating point absorption to magic number
------

\\
First some bookkeeping. I'm going to be assuming double precision (*binary64*) and the default rounding mode (nearest and ties to even). We'll start by considering how *effective* additions and subtractions behave.  Given the examples above:

    s = a + b;
    d = a - b;

We'll call `s` an *effective* addition if `a` and `b` have the same signs and otherwise an *effective* subtraction.
We'll call `d` an *effective* addition if `a` and `b` have opposite signs and otherwise an *effective* subtraction.


So we can consider these two cases individually and under the assumptions that both are positive, ordered and ignore zero: $ a\ge b \gt 0$

Floating point absorption (for addition/subtraction) is when the smaller input doesn't contribute to the result.  In the effective addition case we have `(a+b)==a`. Let's consider a couple of examples using some 4-bit floating point format where on of the left $b$ is the smallest normal in the format and $a$ is $2^p~b$  (where $p$ is the number of precision bits in the format so 4 in this example). On the right we change $b$ to the largest denormal:


    ..1000....     ..1000....     a
    ......1000     ......0111     b
    ..10001000     ..10000111    a+b (before rounding)
	      ^              ^       location of rounding unit
	..1000....     ..1000....    a+b (rounded)
		  

In both of these cases $b$ doesn't contribute (is absorbed). In the left hand case there's no rounding because we have a tie and the lowest bit of $a$ is zero.  If lowest bit of $a$ had been set or if any other bit of $b$ was set then $b$ would have contributed but that's okay since both are normal numbers.  In the right hand case $b$ is absorbed by the example $a$ and all larger numbers.  This is our magic number.

For the *effective* subtraction case we need to borrow from the larger so let's look at two examples with $b$ set to the max denormal and $a$ as above and one ULP larger.

    ..1000....     ..1001....     a
    ......0111     ......0111     b
    ...1111001     ..10001001    a-b (before rounding)
	       ^             ^       location of rounding unit
	...1111...     ..1001....    a-b (rounded)
	
In the left hand case (which was good above) $b$ does contribute. However if $a$ is the next representable FP number then it is absorbed.

So we can describe the magic number as:  if $\text{min}$ is the smallest normal number of a $p$ precision binary floating point number format then the threshold $T$ is:

$$
T = 2^p~\text{min}
$$

In binary64 we have $min=2^{-1022}$ and $p=2^{53}$ so $T=2^{-969}$

For the effective addition case the test value must be greater-than or equal to $T$ and the effective subtraction case (or don't know which) it must be greater than. 

Recap of the `threshold` value and test:
1. $\text{min}$ allows enough space for the smaller magnitude input to be the any number with the smallest representable exponent. If it's a denormal then what would be highest bit of the previous is zero.
2. multiply by $2^p$ to allow enough space such that the larger magnitude input's bits can't overlap with those of a number with the smallest representable exponent.
3. a basic version of the test becomes:
   * if can only be an effective addition then: `if (|r| >= threshold) { expected common case }`
   * otherwise `if (|t| >  threshold) { expected common case }`
4. the above test doesn't capture all $a$ and $b$ considered legal. only those that can be determined by examining the result of the add/sub.
5. as noted at the start the test *can* be used to catch special inputs and overflows

<br>

------

Some extensions and bells and whistles
------

\\
Before stopping to type I'll toss out a few quick extensions of the base idea.

First let's figure out the magic number for double-doubles. Briefly: each value is represented by a hi/lo pair such that (hi+lo == hi). So:
1. for the smaller magnitude we need the term $2^{53}~\text{min}$ so there's enough space for both words to be normal numbers.
2. for the larger magnitude we need to allow space for it to not overlap giving an additional $2^{106}$ factor
3. giving:  $T = 2^{106} \cdot 2^{53} \cdot 2^{-1022} = 2^{-863}$


Next let's revert back to plain double computations and pretend the next operation (after the add or sub) is a square root and the sum/diff can be negative. Well we can obviously drop the `fabs` on the add/sum since a negative result doesn't need to be directed to the common case but let's pretend we have a good reason to want to include the root in the optimistic computation prior to the comparison. Well the square root is strictly increasing so if $r > T$ then $\sqrt{r} > \sqrt{T}$ over Reals.  Over floating point we might need to tweak $\sqrt{T}$ to account for how our magic cutoff rounds. On average square roots send two inputs to a given output map since the interval $\left[1,4\right)$ maps to $\left[1,2\right)$. Let's peek at the correctly rounded square of $T$ and a couple of its successors:

    RN(sqrt(0x1.0000000000000p-969)) = 0x1.6a09e667f3bcdp-485  // sqrt(T)
    RN(sqrt(0x1.0000000000001p-969)) = 0x1.6a09e667f3bcdp-485  // sqrt(T+ulp)
    RN(sqrt(0x1.0000000000002p-969)) = 0x1.6a09e667f3bcep-485  // sqrt(T+2 ulp)

(FWIW: the input power is odd so sqrt(2) and exp adjust). So it happens that the correctly rounded square of $T$ is fine for a greater than comparison.

	

