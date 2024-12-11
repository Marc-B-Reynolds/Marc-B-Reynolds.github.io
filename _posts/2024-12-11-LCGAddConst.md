---
layout:       post
title:        Comments on parameterized SplitMix and PCG style generators
categories:   [math]
tags:         [random]
description:  'Builds up a rational for a presented method for generating the additive constant of a LCG or Weyl sequence based PRNG.'
plotly:       false
---

This is minimal summary of a blog post I've never pushed myself to complete on parameterized *SplitMix* and [PCG](https://en.wikipedia.org/wiki/Permuted_congruential_generator) type PRNGs. 
It's half rambling because I'm quickly tossing it out on a whim.

Roughly:

1. Parameterization is sound provided testing has been performed on the *fundamental* sequence (see below).
2. Provided the first step has been performed then any reasonable choice of additive constant cannot be weaker.

This is stripped down to give minimal justification for a proposed method to generate the additive constant given at the end. The rest is just justicification of why I think that method and parameterization are sound.

<br>

------

State update property sketch
------

\\
First let's note how the states of both are updated:

1. *SplitMix* walks an additive recurrence (integer Weyl sequence)
2. *PCG* walks a linear congruential generator (LCG)

then both apply a *bit finalizer* to return a result.

Both of these state sequences (Weyl & LCG) have a *fundamental* sequence:

1. Weyl is the counting numbers: {0,1,2,3,...}
2. LCG is the chosen multiplicative constant (`M`) with the additive constant set to '1'.

Calling the chosen additive constant of either `A` then transforming the fundamental sequence into the one actually used is simply multipling each produced element by `A`. So we can consider both to be walking their fundamental sequence and the integer multiply by `A` to be a "for free" integer bijection step.  So reasonable choices of $A$ are those we'd find used in finalizer steps whose most basic properties are those of the average random integer (popcount approximately half the bit width being the most basic).

Repeating the above using math for the LCG case if we call the fundamental sequence $y$ and a general sequence $x$ both with the first elements $y_0$ and $x_0$ set to zero then we have the following recurrence relations for both and finally the closed-form transform of $y_i$ to $x_i$:

$$ 
\begin{eqnarray*}
\mathbf{y}_0     & = & 0                     \\
\mathbf{x}_0     & = & 0                     \\
\\
\mathbf{y}_{i+1} & = & M~\mathbf{y}_{i} + 1 \\
\\
\mathbf{x}_{i+1} & = & M~\mathbf{x}_{i} + A \\
\mathbf{x}_i     & = & A~\mathbf{y}_{i}
\end{eqnarray*}
$$

<details markdown="1">
<summary>(click to expand) code snippets just to be certain this is clear:</summary>

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ c
#include <stdint.h>

// generic state data & parameterized constant for both
typedef struct { uint64_t state; uint64_t inc; } prng_t;

// LCG multiplicative constant
static const uint64_t prng_mul_k = UINT64_C(0xd1342543de82ef95);

// example addtive constant for both: scaled rounded to odd 1/phi
// only hard requirement is an odd number
static const uint64_t prng_add_k = UINT64_C(0x9e3779b97f4a7c15);

// example "standard" parameterized state updates

static inline uint64_t weyl_update(prng_t* rng)
{
  uint64_t state = rng->state;
  
  rng->state += rng->inc;

  return state;
}

static inline uint64_t lcg_update(prng_t* rng)
{
  uint64_t state = rng->state;
  
  rng->state = prng_mul_k * state + rng->inc;

  return state;
}

// equivalent parameterized state updates

// always walking this instead of the above:
// 1) 'state' walks the fundamental sequence (inc=1)
// 2) but it returns the same result

static inline uint64_t weyl_update_equivalent(prng_t* rng)
{
  uint64_t state = rng->state;
  
  rng->state += 1;

  // logically a bit mixing/finalizing/hashing step
  return rng->inc * state;
}


static inline uint64_t lcg_update_equivalent(prng_t* rng)
{
  uint64_t state = rng->state;
  
  rng->state = prng_mul_k * state + 1;

  // logically a bit mixing/finalizing/hashing step
  return rng->inc * state;
}

// mini sanity check

#include <stdio.h>

#define TRIALS 0xfffff              // doesn't matter

int main(void)
{
  prng_t lcg0  = {.state=0,.inc=prng_add_k};
  prng_t lcg1  = {.state=0,.inc=prng_add_k};
  prng_t weyl0 = {.state=0,.inc=prng_add_k};
  prng_t weyl1 = {.state=0,.inc=prng_add_k};

  for(uint32_t i=0; i<TRIALS; i++) {
    uint64_t l0 = lcg_update(&lcg0);
    uint64_t l1 = lcg_update_equivalent(&lcg1);
    uint64_t w0 = weyl_update(&weyl0);
    uint64_t w1 = weyl_update_equivalent(&weyl1);

    if ((l0 == l1) && (w0 == w1)) continue;

    printf("never reached\n");
    return -1;
  }
  return 0;
}

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

</details>

<br>

------

Sketch of implications of the above properties
------

\\
This section is to minimally expand on the previous with some extra details not really needed for the proposed code and spitball a few side topics.

Starting with *SplitMix64* we have an effective state update of simply walking counting numbers, the implied integer product and the application of a bit finalizer passing statistical randomness testing. An example of independently performed testing by Daniel Lemire can be found in this [blog post](https://lemire.me/blog/2017/08/22/testing-non-cryptographic-random-number-generators-my-results/). Counting numbers are the most unrandom sequence possible.

PCGs up the ante by effectively walking some fundamental linear congruential sequence (I babble about some properties [here](https://marc-b-reynolds.github.io/math/2017/09/12/LCS.html)), the implied product and finalizer.  In some respects power-of-two LCGs have a worse reputation than they deserve. The have two main problems. The first in "in the planes" which I'm going to ignore because one should choose the multiplicative constant from a well chosen [table](https://gist.github.com/Marc-B-Reynolds/53ecc4f9648d9ba4348403b9afd06804) which are designed to minimize this problem. The second being the bit period properties which make low bit useless. As an example some 128-bit LCGs that only return the high 64-bits have been shown to pass the *BigCrush* suite of *TestU01*.

The addition of the bit finalizer transforms their weaknesses into strengths. In both cases we're transforming simple linear sequence with very well known properties. Notably how both behave on some power-of-two interval which in combination with using a well performing bit finalizer at full input to output bit [strict avalanche criterion](https://marc-b-reynolds.github.io/math/2019/08/10/Avalanche.html) testing gives a reasonable expection of how all same sized power-of-two intervals behave across the whole sequence. Not a replacement for statistical testing but additional confidence. George Marsagila called all LCGs with same multiplicative constant "statistically equivalent" because of their relation to the fundamental sequence and here we have all power-of-two subsequences are statistically equivalent under the assumption of full avalanche.

A few more asides. The properties suggests a component-by-component by testing strategy for designing a PCG variant: 
1. finalizer tested in isolation
2. finalizer mixing a Weyl sequence 
3. finalizer mixing fundamental LCG

I'm mentioning this because since *SplitMix64* passes statistical testing I use its finalizer ([MIX14](http://zimbry.blogspot.com/2011/09/better-bit-mixing-improving-on.html)) for PCG variants.

<details markdown="1">
<summary>(click to expand) example code snippet</summary>

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ c

static const uint64_t prng_mul_k = UINT64_C(0xd1342543de82ef95);
static const uint64_t prng_add_k = UINT64_C(0x2545f4914f6cdd1d);

typedef struct { uint64_t state; uint64_t inc; } prng_t;

static inline uint64_t prng_mix_64(uint64_t x)
{
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
  x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
  x = (x ^ (x >> 31));
  return x;
}

// not mentioned in text: since when returning 32-bit integers,
// single precision floats, etc. we're going to discard at least
// 32-bits so with a finalizer whose last step only effect
// the bottom 32-bits..well that step can be dropped provided
// they are discarded. The downside of this is that we have
// another "sequence" that needs to be run through PRNG testing.
// This example has "not" so it's only for illustrative purposes.
// But it's fine for causal random number generation as it's as 
// strong a finalizer as MIX14. 
static inline uint64_t prng_mix_32(uint64_t x)
{
  // xxhash finalizer
  x = (x ^ (x >> 33)) * 0xc2b2ae3d27d4eb4f;
  x = (x ^ (x >> 29)) * 0x165667b19e3779f9;
//x = (x ^ (x >> 32));                         // dropped step

  return x;  // caller must discard at least 32 low bits
}

static inline uint64_t prng_u64(prng_t* prng)
{
  uint64_t s  = prng->state;
  uint64_t r  = prng_mix_64(s);

  prng->state = prng_mul_k * s + prng->inc;
  
  return r;
}

static inline uint32_t prng_u32(prng_t* prng)
{
  uint64_t s  = prng->state;
  uint32_t r  = (uint32_t)(prng_mix_32(s)>>32);

  prng->state = prng_mul_k * s + prng->inc;
  
  return r;
}

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
</details>

<br>

------

A parameterized additive constant generation scheme
------

\\
Here's the proposed additive constant generation scheme. The goal is generate likely "candidates" for a integer multiplicative constant as cheaply as possible in a way that is concurrent friendly and easy to use.

1. Walk a simple counter (`c`):
   * The initial value of the counter to be set to anything for simplicity which would produce generators with distinct sequences per run. Say if loaded and stored per run then we'll never-ever see the same *generators* twice (see below). Note this is independent from the common technique of randomizing the state which changes our index in the sequence.
   * convert `c` to odd (note this discards the top bit as written).
   * multiply `c` by a *magic constant* to convert it into the walking the odd elements of a Weyl sequence to produce a candidate value. A Weyl sequence because it's cheap and as a (weak) [low-discrepancy sequence](https://en.wikipedia.org/wiki/Low-discrepancy_sequence) will produce values that (roughly) are *uniformly* scattered across the set of integers.
2. Reject values that have a popcount outside some window of half the number of bits (goto 1 on fail)
3. Reject values that have too few runs of 1s with respect to the population count as a cheap approximation to the expected number & lengths of bit runs. We're only trying to avoid obviously bad choices. (goto 1 on fail)
4. Accept and return the candidate.

The number of 64-bit integers with popcounts within window $\left(w\right)$ of 32 is:

$$ 
\text{total}\left(w\right) = \binom{64}{32} + 2 \sum_{n=1}^{w} \binom{64}{32-i}
$$

with an unreasonable window of 0 this is approximately $2^{60.67}$ distinct generators and with the example choice in the code below it's $\approx\left(2^{63.95}\right)$. The next filtering (low number of 1 run counts) rejects approaching zero of the previous so doesn't effect these approximate numbers. 

Example implementation:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ c
// acceptable window for population count:
//  number of 64-bit integers with popcount within 'w' of 32:
//    total[w] = Binomial[64,32] + 2 Sum[Binomial[64,32-i], {i,1,w}]
//  expected number of iterations to find integer within 'w':
//   2^64 / total[w]
//  seen : average iterations seen in 2^30 trials
//  max  : max     iterations seen in 2^30 trials
//
//  w   expected  seen  max
//  1 : 3.42443 : 3.37 : 61
//  2 : 2.13816 : 2.11 : 32
//  3 : 1.61742 : 1.60 : 21
//  4 : 1.35215 : 1.33 : 16
//  5 : 1.20285 : 1.19 : 12
//  6 : 1.11535 : 1.11 : 10
//  7 : 1.06376 : 1.06 :  8
//  8 : 1.03388 : 1.03 :  7
//  9 : 1.01715 : 1.02 :  7
// 10 : 1.00821 : 1.01 :  5

// compile time choice of acceptable window of
// population count
static const uint64_t prng_internal_inc_pw = 8;

// MAGIC! (I made it up). 
// SRLY: not off-the-shelf so unlikely to be used anywhere else
// ditto for not being a 64-bit prime with high bit set
// square-free: 3*5*11*601*114967123673219
// reasonable bit distribution (pop=38):
// 1001111000110111011110011011100101101111010010100111100010010111
// 
// Off the shelf reasonable choices include the classic scaled golden
// ratio and similar constants that behave like a weak LDS. 
static const uint64_t prng_internal_inc_k = 0x9e3779b96f4a7897;

// simple counter. can be init time set to anything
// so that generators sequences aren't the same each run.
// (even zero as the routine well be rejected it). top bit is
// discarded in converting to odd.
static _Atomic uint64_t prng_internal_inc_id = 1;

uint64_t lcg_additive_next(void)
{
  uint64_t b;

  do {
    // atomically increment the global counter and
    // convert it into a candidate additive constant
    b  = atomic_fetch_add_explicit(&prng_internal_inc_id,
				   1,
				   memory_order_relaxed);
    b  = (b<<1)|1;            // make odd
    b *= prng_internal_inc_k; // convert to n^th element

    // compute the population count and test value 't'
    uint32_t pop = pop_64(b);
    uint32_t t   = pop - (32-prng_internal_inc_pw);

    // is the popcount within the specified window?
    //   pop in [32-w,32+w]
    if (t <= 2*prng_internal_inc_pw) {

      // don't want to do anything too complicated so
      // just compute the number of runs of 1s in 'b'
      uint32_t str = pop_64(b & (b ^ (b >> 1)));

      // this pretty much always passes. just prevents
      // too many long runs of 1s. requires that runs
      // (on average) are less than 4 bits long.
      if (str >= (pop >> 2)) return b;
    }
  } while(1);
}
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

