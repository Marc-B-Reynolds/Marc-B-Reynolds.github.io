---
layout:       post
title:        Can an involution be a competitive bit finalizer?
categories:   [math]
tags:         [integer,hash,random]
description:  A feasibility study to determine if an involution can be competitive (quality and cost) bit mixer.
---

\\
You know if someone had asked me this question a few months ago I would have responded that it's was very unlikely and my facial expression would have looked like the person was eating huge globs of mayonnaise straight out of the jar.  Let's jump backward and ballpark some defs:

* An [involution](https://en.wikipedia.org/wiki/Involution_(mathematics)) is a function that's its own inverse:  $f\left(f\left(x\right)\right) = f^2\left(x\right) = x$
* A [bijection](https://en.wikipedia.org/wiki/Bijection) is a function that has an inverse:  $f^{-1}\left(f\left(x\right)\right) = f\left(f^{-1}\left(x\right)\right) = x$

where the $n$ in $f^n$ is denoting $n$ [iterations](https://en.wikipedia.org/wiki/Iterated_function) of $f$. And I'm using *"bit finalizer"* to mean a bit mixing function like that typically used as the final step of a hash function. I'm going to be using the terminology and heat map figures from the previous post [*"Comments on the Avalanche Effect"*](http://marc-b-reynolds.github.io/math/2019/08/10/Avalanche.html).

Now back to the beginning. I have a couple of issues with involutions as bit finalizers. The first is there aren't many problems that need such a thing (roughly encoding and decoding are the same thing). The second issue it's more computational work to hit reasonable properties but less that I would have previously thought so here we are...me typing.

Code for this post can be found [**here**](https://github.com/Marc-B-Reynolds/Stand-alone-junk/blob/master/src/Posts/involution_mix.c) and for quick eye-balling of functions on [**godbolt**](https://gcc.godbolt.org/z/f_Sds4).
<br>

------

The set up
------

\\
What got me thinking about this problem was a question on twitter by @paniq where he wanted to pack a 2D coordinate in a integer, perform an involution and unpack back to 2D. This would allow to "statelessly" swap the data between pairs of coordinates for an entire grid in a random manner. My suggested attempt at a solution was something like this:

{% highlight c %}
// bad idea #1
uint32_t g0(uint32_t x)
{
  if (x & 1) {
    // call this f(x)
    x *= 0xac564b05;
    x += 0x85ebca77;
  }
  else {
    // and this f^{-1}(x)
    x -= 0x85ebca77;
    x *= 0xdc33c9cd;
  }
  return x;
}
{% endhighlight %}

\\
This is a pair of bijections (one's the inverse of the other) in the same form as a [LCG](https://en.wikipedia.org/wiki/Linear_congruential_generator). These have the know property that they map odd to even and even to odd. By inspecting the value of the low bit (determine if odd or even) we can choose $f\left(x\right)$ for odd and $f^{-1}\left(x\right)$ for even inputs and this composite function is an involution. The "problem" of @paniq has two additional wrinkles. First we functions that work on the same size of the texture being manipulated and small sized make finding a good bit mixer harder. As an example if the texture were 256x256 we'd need an involution for a 16-bit integer. The second trick is the function needs to be parameterized so this swapping operation can been performed differently in some number of way...all of which need to be somewhat independent. These move the problem out of the territory of "fun little thing to think about" to "that sounds like a lot of work". So I'm going to blissfully ignore these extra requirements.

The visualization of the strict avalanche criterion (**SAC**) bias of the above function:

<p align="center" vertical-align="middle">
<canvas id="g0_le"></canvas>
<canvas id="g0_he"></canvas>
</p>

\\
The left hand is sampled using Sobol and right hand a pseudo-random sequence.

<br>

------

The plot thickens
------

\\
Harold Aptroot spiced things up by noting that a similarity transform of an involution is an involution. Sweet! This is something I can work with. Let's run through the justification in painful details using a matrix like notation where $A$ is a bijection and $B$ an involution:

$$
\begin{align*}
  \left( ABA^{-1}\right)^2
  & = \left( ABA^{-1} \right) \left( ABA^{-1} \right) \\
  & =  AB \left( A^{-1}A \right) BA^{-1} \\
  & =  ABBA^{-1}  \\
  & =  AB^2A^{-1} \\
  & =  AA^{-1} \\
  & =  I
\end{align*}
$$

Bullet point it time:

* We're saying some function $ABA^{-1}$ is an involution which has a period of 2 so squaring it should yield the identity $I$. Start expanding.
* $A$ is a bijection so it has an inverse $A^{-1}$, multipled together gives $I$ which reduces to a 70s Swedish pop group.
* $B$ is an involution, squaring it reduces to $I$ and finally $A$ and it's inverse cancel each other out.

<br>

------

Okay..but does it **do** anything?
------

\\
The rest of this is just a quick feasibility study to see if we can construct an involution that is competitive against standard methods of bit mixing both in terms of statistical measures and runtime costs. It is pretty much a dump of what I walked through to get a feel for the answer. No real search methods were used..just me hacking stuff and seeing what happened.

First we need a fast involution to be transformed and hopefully one that performs some bit mixing as well. The most obvious choice (at least to me) is a right xorshift by at least half the bit-width of the register. Since we're using a bit operation $\left(\mathbb{F}_2\right)$ for the involution I want to choose a modulo integer operation for the bijection (alternating fields between each logical operation). So about the cheapest thing we can do here is multiply by an odd integer (the inverse of which can be found by [mod inverse](http://marc-b-reynolds.github.io/math/2017/09/18/ModInverse.html)). So I simply grabbed an odd constant from L'Ecuyer's [table paper](http://www.ams.org/journals/mcom/1999-68-225/S0025-5718-99-00996-5/S0025-5718-99-00996-5.pdf) and tweaked the xorshift amount until the bias was low.


{% highlight c %}
#define M0 0x5f356495
#define M1 0x32c446bd

uint32_t f0(uint32_t x)
{
  x *= M0;
  x ^= (x>>25);
  x *= M1;
  return x;
}
{% endhighlight %}


<p align="center" vertical-align="middle">
<canvas id="f0_le"></canvas>
<canvas id="f0_he"></canvas>
</p>

\\
This is doing suprisingly well for a short sequence involution. And again the above constants are going to be sub-optimal as are all the functions that follow.

Instead of beefing up the bijection I wanted to see how changing the involution to an [xor-rotate](http://marc-b-reynolds.github.io/math/2017/10/13/XorRotate.html) would play out. First we need some helper functions:

{% highlight c %}
inline uint32_t rot(uint32_t x, uint32_t i)
{
#if defined(__clang__)
  return __builtin_rotateleft32(x,i);
#elif defined(_MSC_VER)
  return _rotl(x,i);
#else
  return (x << i) | (x >> (-i & 31));
#endif
}

// I + C^a + C^b
inline uint32_t xor_rot2(uint32_t x, uint32_t a, uint32_t b)
{
  return x^rot(x,a)^rot(x,b);
}

{% endhighlight %}

\\
I kept the multiplicative constants the same and tweaked the xor-rotate choice:

{% highlight c %}
uint32_t f1(uint32_t x)
{
  x *= M0;
  x = xor_rot2(x,6,22);
  x *= M1;
  return x;
}
{% endhighlight %}

<p align="center" vertical-align="middle">
<canvas id="f1_le"></canvas>
<canvas id="f1_he"></canvas>
</p>

\\
This was a nice little jump, we're still on the cheap side and at this point I'd seen about enough. To kill off the bias on the low bits I just added (didn't bother to tweak) right xorshift by half the bitwidth.

{% highlight c %}
uint32_t f2(uint32_t x)
{
  x ^= x >> 16;
  x *= M0;
  x = xor_rot2(x,6,22);
  x *= M1;
  x ^= x >> 16;
  return x;
}
{% endhighlight %}


<p align="center" vertical-align="middle">
<canvas id="f2_le"></canvas>
<canvas id="f2_he"></canvas>
</p>

\\
At this point we're looking good and we're also at about exactly the same runtime cost as current standard methods that are structured like this:


{% highlight c %}
// 32-bit MurmurHash3, xxHash et al. look like this:
uint32_t common_bijective_finalizers(uint32_t x)
{
  x ^= x >> S0; x *= K0;
  x ^= x >> S1; x *= K1;
  x ^= x >> S2;
  return x;
}
{% endhighlight %}

\\
As a final "I want to see what happens" test I swapped out the outer xorshift with an xor-rotate randomly chosen from those that have a period of 32 and a 3 rotate inverse:

{% highlight c %}

// C^a + C^b + C^c
inline uint32_t xor_rot3(uint32_t x, uint32_t a, uint32_t b, uint32_t c)
{
  return rot(x,a)^rot(x,b)^rot(x,c);
}

uint32_t f3(uint32_t x)
{
  x = xor_rot2(x,11,16); x *= M0;
  x = xor_rot2(x, 6,22); x *= M1;
  x = xor_rot3(x,10,21,26);
  return x;
}
{% endhighlight %}

<p align="center" vertical-align="middle">
<canvas id="f3_le"></canvas>
<canvas id="f3_he"></canvas>
</p>

<br>

------

Crush-ing some numbers
------

\\
This section spews out the results you'd see from running the example code. It measures SAC using three different sampling methods:

1. $\text{cn}$: counting numbers (well programmer version thereof aka starting from zero) for low entropy input. 
2. $\text{ss}$: Sobol sequence to get good coverage of input space.
3. $\text{he}$: random numbers for high entropy input.

\\
Each set is sampled $2^{23}$ times. The **max bias** measured as described in the *Avalanche Effect* post then multipled by 100 to give a percent bias of the worst case. Also a goodness-of-fit (**GOF**) test is performed to get a global estimate. Specifically [Pearson's chi-squared](https://en.wikipedia.org/wiki/Pearson%27s_chi-squared_test) is computed, takes the square root and multiplied by 100.

Additionally three batteries of [TestU01](http://simul.iro.umontreal.ca/testu01/tu01.html) are run on the input. Each hash function is treated like a random number generator where we walk through counting numbers and apply the hash.

1. **Rabbit**: computes 33 statistics from 26 base tests drawing approximately a user defined number of sample (specified in bits). My runs draw 1024 samples of the function so is setup for 1024*32 bits.
2. **SmallCrush**: computes 15 statistics from 10 base tests drawing approximately $51,320,000$ samples
3. **Crush**: computes 144 statistics from 96 base tests drawing approximately $2^{35}$ samples

The reported numbers are how many statistics of the battery failed.

For reference I've include 32-bit finalizer of MurmurHash3, xxHash and a pair by [Chris Wellons](https://github.com/skeeto/hash-prospector). All except *triple32* are two xorshift-multiplies followed by an xorshift and *triple32* add an additional xorshift-multiply. 


Clicking on the headers sorts the data.

<br>

{: .center }
| hash      | max bias (cn) | max bias (ss) | max bias (he) | GOF (cn)| GOF (ss) | GOF (he) | Rabbit| SmallCrush | Crush |
| :---:     |  ---:      |  ---:    |  ---:    |  ---:    |  ---:    | ---:    | ---:  | ---:       |---:       |
| murmur3   |  0.229263  |  0.518417  |  0.207162  |  0.052966  |  0.092238  |  0.043021  |     1 |        2   | 24 |
| xxhash32  |  0.377083  |  0.579166  |  0.433731  |  0.069322  |  0.090209  |  0.066725  |     1 |        2   | 39 |
| triple32  |**0.135088**|**0.156140**|  0.130367  |**0.044136**|**0.045361**|**0.034147**|     1 |        0   |  5 |
| lowbias   |  0.169849  |  0.266051  |**0.122666**|  0.047634  |  0.068301  |  0.040458  |     2 |        2   | 44 |
| f2        |  0.409937  |  0.393176  |  0.398207  |  0.054149  |  0.070380  |  0.045054  |     1 |        0   | 13 |
| f3        |  0.591612  |  0.445747  |  0.516152  |  0.056496  |  0.052190  |  0.045781  |     1 |        0   | **3** |


\\
The three inital functions are failures but here's their data as well (minus Crush since it's a waste of time).

{: .center }
| hash      | max bias (cn) | max bias (ss) | max bias (he) | GOF (cn)| GOF (ss) | GOF (he) | Rabbit| SmallCrush | 
| :---:     |  ---:    |  ---:    |  ---:    |  ---:    |  ---:    | ---:    | ---:  | ---:       |
| g0        |100.000000|100.000000|100.000000| 76.090304|82.050457 |76.090831|     7 |       15   |
| f0        |100.000000|100.000000|100.000000| 23.667056|38.865972 |23.663481|     8 |       14   |
| f1        |100.000000|100.000000|100.000000| 20.454587|20.904367 |20.456410|     1 |        7   |


<br>

------

Conclusions
------

\\
I normally don't have a conclusion section. I talk some trash, you can dig through with a stick and decide if there's anything you want to take home to try out. In this post I goofing around abit so I'm going to toss out my take-a-way.

* The answer to the title is yes.

\\
Bonus takeaways:
* xor-rotates might be worth investigating for various bit mixing operations
* Chris Wellons' finalizers perform rather well

<br>

------

Crush output <small>...so you don't have to</small>
------

\\
Running the Crush battery takes about an hour so to save you the PITA here is what you would get:

```
========= Summary results of Crush =========

 Version:          TestU01 1.2.3
 Generator:        f1
 Number of statistics:  144
 Total CPU time:   00:26:31.17
 The following tests gave p-values outside [0.001, 0.9990]:
 (eps  means a value < 1.0e-300):
 (eps1 means a value < 1.0e-15):

       Test                          p-value
 ----------------------------------------------
  1  SerialOver, t = 2                eps  
  2  SerialOver, t = 4               9.4e-4
  3  CollisionOver, t = 2           6.2e-11
  4  CollisionOver, t = 2             eps  
  6  CollisionOver, t = 4             eps  
  8  CollisionOver, t = 8             eps  
 10  CollisionOver, t = 20            eps  
 11  BirthdaySpacings, t = 2          eps  
 13  BirthdaySpacings, t = 4       2.4e-227
 15  BirthdaySpacings, t = 7          eps  
 16  BirthdaySpacings, t = 8          eps  
 17  BirthdaySpacings, t = 8          eps  
 18  ClosePairs NP, t = 2            6.7e-4
 18  ClosePairs mNP1, t = 2          2.9e-6
 18  ClosePairs mNP2, t = 2          7.3e-5
 18  ClosePairs NJumps, t = 2        2.7e-5
 23  SimpPoker, d = 16               6.2e-6
 24  SimpPoker, d = 16                eps  
 26  SimpPoker, d = 64                eps  
 28  CouponCollector, d = 4           eps  
 30  CouponCollector, d = 16          eps  
 32  Gap, r = 27                      eps  
 33  Gap, r = 0                      6.8e-7
 34  Gap, r = 22                      eps  
 36  Run of U01, r = 15               eps  
 41  MaxOft, t = 5                  1 - eps1
 41  MaxOft AD, t = 5               1.4e-10
 42  MaxOft, t = 10                 1 - eps1
 43  MaxOft, t = 20                 1 - eps1
 44  MaxOft, t = 30                 1 - eps1
 49  AppearanceSpacings, r = 0        eps  
 50  AppearanceSpacings, r = 20       eps  
 52  WeightDistrib, r = 8           1.1e-12
 54  WeightDistrib, r = 24            eps  
 63  GCD, r = 0                       eps  
 64  GCD, r = 10                      eps  
 65  RandomWalk1 H (L = 90)           eps  
 65  RandomWalk1 M (L = 90)           eps  
 66  RandomWalk1 H (L = 90)           eps  
 66  RandomWalk1 M (L = 90)           eps  
 66  RandomWalk1 J (L = 90)           eps  
 66  RandomWalk1 R (L = 90)           eps  
 66  RandomWalk1 C (L = 90)           eps  
 67  RandomWalk1 H (L = 1000)         eps  
 67  RandomWalk1 M (L = 1000)         eps  
 67  RandomWalk1 R (L = 1000)         eps  
 67  RandomWalk1 C (L = 1000)         eps  
 68  RandomWalk1 H (L = 1000)         eps  
 68  RandomWalk1 M (L = 1000)         eps  
 68  RandomWalk1 J (L = 1000)       4.0e-12
 68  RandomWalk1 R (L = 1000)         eps  
 68  RandomWalk1 C (L = 1000)         eps  
 69  RandomWalk1 H (L = 10000)        eps  
 69  RandomWalk1 M (L = 10000)        eps  
 69  RandomWalk1 R (L = 10000)      1.1e-16
 69  RandomWalk1 C (L = 10000)        eps  
 70  RandomWalk1 H (L = 10000)        eps  
 70  RandomWalk1 M (L = 10000)        eps  
 70  RandomWalk1 R (L = 10000)        eps  
 70  RandomWalk1 C (L = 10000)        eps  
 72  LinearComp, r = 29             1 - 4.7e-14
 74  Fourier3, r = 0                  eps  
 75  Fourier3, r = 20                 eps  
 77  LongestHeadRun, r = 20          7.2e-4
 80  HammingWeight2, r = 0          1 - eps1
 81  HammingWeight2, r = 20         1 - eps1
 82  HammingCorr, L = 30            1 - eps1
 83  HammingCorr, L = 300           1 - eps1
 84  HammingCorr, L = 1200          1 - eps1
 85  HammingIndep, L = 30             eps  
 86  HammingIndep, L = 30             eps  
 87  HammingIndep, L = 300            eps  
 88  HammingIndep, L = 300            eps  
 89  HammingIndep, L = 1200           eps  
 90  HammingIndep, L = 1200           eps  
 92  Run of bits, r = 20              eps  
 95  AutoCor, d = 30                  eps  
 96  AutoCor, d = 10                  eps  
 ----------------------------------------------
 All other tests were passed


========= Summary results of Crush =========

 Version:          TestU01 1.2.3
 Generator:        f2
 Number of statistics:  144
 Total CPU time:   00:26:56.26
 The following tests gave p-values outside [0.001, 0.9990]:
 (eps  means a value < 1.0e-300):
 (eps1 means a value < 1.0e-15):

       Test                          p-value
 ----------------------------------------------
  1  SerialOver, t = 2                eps  
 11  BirthdaySpacings, t = 2        6.2e-82
 25  SimpPoker, d = 64               6.1e-8
 41  MaxOft, t = 5                  1 - 5.9e-10
 41  MaxOft AD, t = 5                2.9e-9
 42  MaxOft, t = 10                 1 - eps1
 43  MaxOft, t = 20                 1 - 1.7e-11
 44  MaxOft, t = 30                 1 - eps1
 49  AppearanceSpacings, r = 0      1 - eps1
 52  WeightDistrib, r = 8            2.5e-7
 82  HammingCorr, L = 30            1 -  8.4e-6
 85  HammingIndep, L = 30            5.2e-4
 95  AutoCor, d = 30                9.5e-33
 ----------------------------------------------
 All other tests were passed


========= Summary results of Crush =========

 Version:          TestU01 1.2.3
 Generator:        f3
 Number of statistics:  144
 Total CPU time:   00:27:52.02
 The following tests gave p-values outside [0.001, 0.9990]:
 (eps  means a value < 1.0e-300):
 (eps1 means a value < 1.0e-15):

       Test                          p-value
 ----------------------------------------------
 42  MaxOft, t = 10                 1 - 9.9e-12
 43  MaxOft, t = 20                 1 -  4.6e-6
 44  MaxOft, t = 30                 1 - 2.6e-14
 ----------------------------------------------
 All other tests were passed


========= Summary results of Crush =========

 Version:          TestU01 1.2.3
 Generator:        xxhash32
 Number of statistics:  144
 Total CPU time:   00:27:14.14
 The following tests gave p-values outside [0.001, 0.9990]:
 (eps  means a value < 1.0e-300):
 (eps1 means a value < 1.0e-15):

       Test                          p-value
 ----------------------------------------------
  1  SerialOver, t = 2                eps  
  2  SerialOver, t = 4              8.9e-14
  3  CollisionOver, t = 2          4.2e-281
 11  BirthdaySpacings, t = 2        1.1e-93
 13  BirthdaySpacings, t = 4          eps  
 18  ClosePairs NP, t = 2           2.1e-10
 18  ClosePairs mNP, t = 2          1.0e-60
 18  ClosePairs mNP1, t = 2         1.5e-66
 18  ClosePairs mNP2, t = 2         3.4e-20
 18  ClosePairs NJumps, t = 2       6.5e-87
 23  SimpPoker, d = 16                eps  
 24  SimpPoker, d = 16                eps  
 25  SimpPoker, d = 64                eps  
 26  SimpPoker, d = 64                eps  
 28  CouponCollector, d = 4           eps  
 30  CouponCollector, d = 16         3.7e-9
 31  Gap, r = 0                       eps  
 32  Gap, r = 27                      eps  
 33  Gap, r = 0                     1.1e-14
 34  Gap, r = 22                      eps  
 35  Run of U01, r = 0              6.8e-14
 36  Run of U01, r = 15               eps  
 41  MaxOft, t = 5                  1 - eps1
 42  MaxOft, t = 10                 1 - eps1
 43  MaxOft, t = 20                 1 -  1.4e-9
 44  MaxOft, t = 30                 1 - 3.6e-15
 48  SampleCorr                      6.8e-7
 49  AppearanceSpacings, r = 0        eps  
 52  WeightDistrib, r = 8            1.9e-4
 54  WeightDistrib, r = 24          1.3e-12
 57  MatrixRank, 60 x 60             0.9993 
 81  HammingWeight2, r = 20         1 -  3.7e-6
 82  HammingCorr, L = 30            1 -  4.4e-5
 84  HammingCorr, L = 1200          1 -  7.2e-6
 85  HammingIndep, L = 30           1.5e-14
 92  Run of bits, r = 20             8.3e-4
 95  AutoCor, d = 30                1.6e-31
 96  AutoCor, d = 10                2.8e-25
 ----------------------------------------------
 All other tests were passed


========= Summary results of Crush =========

 Version:          TestU01 1.2.3
 Generator:        murmurhash3
 Number of statistics:  144
 Total CPU time:   00:27:08.43
 The following tests gave p-values outside [0.001, 0.9990]:
 (eps  means a value < 1.0e-300):
 (eps1 means a value < 1.0e-15):

       Test                          p-value
 ----------------------------------------------
  1  SerialOver, t = 2                eps  
  3  CollisionOver, t = 2             eps  
  4  CollisionOver, t = 2            6.7e-7
 11  BirthdaySpacings, t = 2        4.7e-42
 13  BirthdaySpacings, t = 4          eps  
 18  ClosePairs NP, t = 2           1.4e-15
 18  ClosePairs mNP, t = 2         3.2e-157
 18  ClosePairs mNP1, t = 2        8.9e-284
 18  ClosePairs mNP2, t = 2        8.0e-123
 18  ClosePairs NJumps, t = 2         eps  
 21  ClosePairsBitMatch, t = 2       4.4e-4
 25  SimpPoker, d = 64                eps  
 26  SimpPoker, d = 64                eps  
 31  Gap, r = 0                       eps  
 33  Gap, r = 0                     2.8e-11
 34  Gap, r = 22                      eps  
 36  Run of U01, r = 15              1.4e-7
 41  MaxOft, t = 5                  2.8e-37
 42  MaxOft, t = 10                 1 - eps1
 43  MaxOft, t = 20                 1 - 1.2e-11
 44  MaxOft, t = 30                 1 - 1.1e-14
 49  AppearanceSpacings, r = 0      1 - eps1
 53  WeightDistrib, r = 16           1.8e-5
 54  WeightDistrib, r = 24            eps  
 ----------------------------------------------
 All other tests were passed


========= Summary results of Crush =========

 Version:          TestU01 1.2.3
 Generator:        lowbias
 Number of statistics:  144
 Total CPU time:   00:27:03.56
 The following tests gave p-values outside [0.001, 0.9990]:
 (eps  means a value < 1.0e-300):
 (eps1 means a value < 1.0e-15):

       Test                          p-value
 ----------------------------------------------
  1  SerialOver, t = 2                eps  
  2  SerialOver, t = 4                eps  
  3  CollisionOver, t = 2             eps  
  4  CollisionOver, t = 2           7.4e-33
 11  BirthdaySpacings, t = 2       4.6e-170
 13  BirthdaySpacings, t = 4          eps  
 15  BirthdaySpacings, t = 7       1.6e-156
 17  BirthdaySpacings, t = 8        2.3e-53
 18  ClosePairs NP, t = 2           1.2e-13
 18  ClosePairs mNP, t = 2         3.2e-157
 18  ClosePairs mNP1, t = 2        1.6e-180
 18  ClosePairs mNP2, t = 2         1.3e-70
 18  ClosePairs NJumps, t = 2      3.9e-249
 21  ClosePairsBitMatch, t = 2       1.1e-4
 23  SimpPoker, d = 16                eps  
 24  SimpPoker, d = 16                eps  
 25  SimpPoker, d = 64                eps  
 26  SimpPoker, d = 64                eps  
 27  CouponCollector, d = 4           eps  
 29  CouponCollector, d = 16          eps  
 30  CouponCollector, d = 16          eps  
 31  Gap, r = 0                       eps  
 32  Gap, r = 27                      eps  
 33  Gap, r = 0                       eps  
 34  Gap, r = 22                      eps  
 35  Run of U01, r = 0                eps  
 36  Run of U01, r = 15              4.3e-4
 37  Permutation, r = 0              1.1e-6
 41  MaxOft AD, t = 5               5.6e-11
 42  MaxOft AD, t = 10              6.4e-15
 43  MaxOft AD, t = 20              1 -  4.8e-8
 44  MaxOft, t = 30                 1 - 3.1e-15
 44  MaxOft AD, t = 30              1 -  6.0e-9
 46  SampleProd, t = 30              0.9992 
 47  SampleMean                      4.5e-6
 48  SampleCorr                     1 -  6.5e-9
 49  AppearanceSpacings, r = 0        eps  
 52  WeightDistrib, r = 8             eps  
 55  SumCollector                     eps  
 75  Fourier3, r = 20                5.7e-5
 81  HammingWeight2, r = 20          0.9992 
 84  HammingCorr, L = 1200          1 -  9.2e-5
 85  HammingIndep, L = 30            1.2e-6
 96  AutoCor, d = 10                 1.1e-8
 ----------------------------------------------
 All other tests were passed


========= Summary results of Crush =========

 Version:          TestU01 1.2.3
 Generator:        triple32
 Number of statistics:  144
 Total CPU time:   00:27:30.04
 The following tests gave p-values outside [0.001, 0.9990]:
 (eps  means a value < 1.0e-300):
 (eps1 means a value < 1.0e-15):

       Test                          p-value
 ----------------------------------------------
 41  MaxOft, t = 5                  1 - 1.8e-13
 42  MaxOft, t = 10                 1 -  1.2e-9
 43  MaxOft, t = 20                 1 -  3.2e-9
 44  MaxOft, t = 30                 1 - eps1
 70  RandomWalk1 H (L = 10000)       7.6e-4
 ----------------------------------------------
 All other tests were passed

```



<script>


const cell_w = 5;  // drawn cell size (pixels) w/o padding

function viz_fill(pix,x,y,w, r,g,b)
{
  var sx = (cell_w+1);
  var i  = 2 + sx*(w+x+y+sx*w*y);
 
  i *= 4;

  for(var k=0; k<cell_w; k++) {
    for(var j=0; j<cell_w;j++) {
      pix[i+0] = r;
  	  pix[i+1] = g;
      pix[i+2] = b;
      pix[i+3] = 255;
	  i += 4;
    }
    i += 4*(2+sx*w-sx);
  }
}

function color(v)
{
  var r,g,b;

  if (v == 0) { r=g=b=0; } else
  if (v >= 0) {
     r = (v>>1) + 128; 
     g = v>>2;
     b = v>>1;
   } else {
     v = -v;
     b = (v>>1) + 128; 
     r = v>>2;
     g = v>>1;
 }

  return [r,g,b];
}

function viz_bias(name,data,w,h)
{
  var cw  = (cell_w+1)*w+1;
  var ch  = (cell_w+1)*h+1;
  var id  = document.getElementById(name);
  var ctx = id.getContext('2d');
  var img = ctx.createImageData(cw,ch);
  var pix = img.data;
  var si  = 0;
  
  id.width  = cw;
  id.height = ch;

  pix.fill(128);

  for (var y=0; y<h; y++) {
    for (var x=0; x<w; x++)  {
	  var c = color(data[si++]);
	  viz_fill(pix, x, y, w, c[0],c[1],c[2]);
    }
  }

  ctx.putImageData(img,0,0)
}

// temp hack
function viz_legend(name,w,h)
{
  var y;
  w = (cell_w+1)*w+1;
  h = (cell_w+1)*h+1;

  var s = 2.0/h;

  var id  = document.getElementById(name);

  if (id != null) {
    var ctx = id.getContext('2d');
    var img = ctx.createImageData(w,h);
    var pix = img.data;
	
	id.width  = 2*w;
    id.height = h;
    pix.fill(255);

    for (var y=0; y<=h; y++) {
      var c = color(Math.round(255.5*(1-s*y)));
	  var i = 4*y*w;
	  
      for (var x=0; x<w; x++)  {
        pix[i+0] = c[0];
  	    pix[i+1] = c[1];
        pix[i+2] = c[2];
        pix[i+3] = 255;
	    i += 4;
      }
    }

    ctx.putImageData(img,0,0)
    return img;
  }
  else console.log('opps');
}

// sigh
function viz_dup(name,data)
{
  var id  = document.getElementById(name);
  if (id != null) {
    id.getContext('2d').putImageData(data,0,0);
  }
}


const g0_le = [-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,128,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,127,-64,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,128,-64,159,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,127,-64,159,-48,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,127,-64,159,-48,151,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,127,-64,159,-48,151,-52,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,128,-64,159,-48,151,-52,-153,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,128,-64,159,-48,151,-52,-153,204,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,128,-64,159,-48,151,-52,-153,204,-25,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,127,-64,159,-48,151,-52,-153,204,-25,140,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,127,-64,159,-48,151,-52,-153,204,-25,140,57,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,-1,128,-64,159,-48,151,-52,-153,204,-25,140,57,-99,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,-2,128,-63,159,-48,151,-52,-153,204,-25,140,57,-99,177,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,-36,145,-55,155,-50,153,-51,-153,204,-25,140,57,-99,177,-39,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,-24,139,-58,156,-49,152,-51,-153,204,-25,140,57,-99,177,-39,-147,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,16,120,-68,161,-47,151,-52,-154,204,-25,140,57,-99,177,-39,-147,-201,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,32,112,-72,163,-46,150,-52,-154,204,-25,140,57,-99,177,-39,-147,-201,228,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,-64,159,-48,151,-52,153,-51,-153,204,-26,140,57,-99,177,-39,-147,-201,228,-13,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,127,64,-96,175,-40,147,-54,-154,205,-25,140,57,-99,177,-39,-147,-201,228,-13,134,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,128,-64,159,-48,151,-52,-153,204,-25,140,57,-99,177,-39,-147,-201,228,-13,134,-60,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,255,0,-128,191,-32,143,-56,-155,205,-25,140,58,-99,177,-39,-147,-201,228,-13,134,-60,158,-255,-255,-255,-255,-255,-255,-255,-255,255,255,0,-127,191,-32,143,-56,-155,205,-25,140,58,-99,177,-39,-147,-201,228,-13,134,-60,158,-49,-255,-255,-255,-255,-255,-255,-255,255,255,-255,-255,255,0,128,-64,-159,207,-24,139,58,-99,177,-39,-147,-201,228,-13,134,-60,158,-49,-152,-255,-255,-255,-255,-255,-255,255,-255,255,-255,255,0,128,-64,-159,207,-24,139,58,-99,177,-39,-147,-201,228,-13,134,-60,158,-49,-152,-203,-255,-255,-255,-255,-255,255,255,255,-255,255,0,128,-64,-159,207,-24,139,58,-99,177,-39,-147,-201,228,-13,134,-60,158,-49,-152,-203,229,-255,-255,-255,-255,255,255,-255,-255,255,255,0,-128,-191,223,-16,135,60,-98,176,-39,-147,-201,228,-13,134,-60,158,-49,-152,-203,229,-13,-255,-255,-255,255,-255,255,-255,255,255,-255,-255,-255,255,0,128,64,-96,175,-40,-147,-201,228,-13,134,-60,158,-49,-152,-203,229,-13,-134,-255,-255,255,255,255,-255,255,-255,255,-255,-255,255,0,128,64,-96,175,-40,-147,-201,228,-13,134,-60,158,-49,-152,-203,229,-13,-134,194,-255,255,-255,255,255,255,255,255,-255,-255,255,0,128,64,-96,175,-40,-147,-201,228,-13,134,-60,158,-49,-152,-203,229,-13,-134,194,30,255,-255,-255,255,-255,-255,255,255,-255,255,0,127,64,-32,-16,-8,12,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0];

const g0_he = [-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,128,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,128,0,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,127,0,32,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,127,0,32,-112,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,128,0,32,-112,-32,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,128,0,32,-112,-32,-143,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,127,0,32,-112,-32,-143,46,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,127,0,32,-112,-32,-143,46,100,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,127,0,32,-112,-32,-143,46,100,52,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,128,0,32,-112,-32,-144,46,100,52,101,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,128,0,32,-112,-32,-143,46,100,52,101,77,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,127,0,32,-112,-32,-143,46,100,52,101,77,-89,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,127,0,32,-111,-32,-144,46,100,52,101,77,-89,172,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,128,0,32,-112,-32,-144,46,100,52,101,77,-89,172,2,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,128,0,32,-111,-32,-143,46,100,52,101,77,-89,172,3,-126,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,128,0,32,-111,-32,-144,46,99,52,101,77,-89,172,3,-126,-10,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,128,0,32,-112,-32,-143,46,100,52,101,77,-89,172,2,-126,-10,95,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,128,0,32,-112,-32,-143,46,100,52,101,77,-89,172,2,-126,-10,95,66,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,-1,128,0,32,-112,-32,-143,46,100,52,101,77,-89,172,2,-126,-10,95,66,94,-255,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,0,127,0,32,-112,-32,-143,46,99,52,101,77,-89,172,2,-126,-10,95,66,94,20,-255,-255,-255,-255,-255,-255,-255,-255,-255,255,4,125,0,32,-111,-32,-143,46,100,52,101,77,-89,172,2,-126,-10,95,66,94,20,118,-255,-255,-255,-255,-255,-255,-255,-255,255,8,123,0,33,-111,-32,-143,46,100,52,101,77,-89,172,2,-126,-10,95,66,95,20,118,20,-255,-255,-255,-255,-255,-255,-255,255,48,104,0,38,-109,-32,-143,46,100,52,101,77,-89,172,3,-126,-10,95,66,94,20,118,20,-118,-255,-255,-255,-255,-255,-255,255,-32,144,0,28,-114,-32,-143,46,99,52,101,77,-89,172,2,-126,-10,95,66,94,20,118,20,-117,-186,-255,-255,-255,-255,-255,255,0,128,0,32,-111,-32,-143,46,100,52,101,77,-89,172,2,-126,-10,95,66,94,20,117,20,-118,-186,221,-255,-255,-255,-255,255,0,127,0,32,-112,-32,-143,46,100,52,102,77,-89,172,3,-126,-10,95,66,94,20,117,20,-118,-186,221,-17,-255,-255,-255,255,0,128,0,32,-112,-32,-143,46,100,52,101,77,-89,172,2,-126,-10,95,66,94,20,117,20,-117,-186,221,-17,2,-255,-255,255,255,0,0,64,-96,-32,-144,44,100,53,101,77,-89,172,2,-127,-10,95,66,94,20,118,20,-117,-186,221,-17,2,126,-255,255,-255,255,0,0,-127,-32,-143,48,100,52,102,77,-89,172,2,-126,-10,95,66,94,20,117,20,-117,-186,221,-17,2,126,64,255,-255,-255,255,0,0,32,-16,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

const f0_le = [-255,-255,-255,-255,-255,-255,255,-16,135,-56,-155,-205,230,12,121,67,-94,175,-40,-148,201,27,-114,-185,220,18,-119,-187,221,-17,136,-59,-255,-255,-255,-255,-255,255,0,56,0,-104,-179,217,0,54,0,-14,134,0,-94,174,0,-44,-149,202,0,-51,-153,204,0,77,0,108,-255,-255,-255,-255,255,0,-127,48,104,-72,163,0,-128,-59,127,-47,0,127,54,0,127,9,123,0,127,1,127,0,-127,32,-113,1,-255,-255,-255,255,0,-127,-191,88,64,80,0,-127,191,-23,-64,-117,128,-64,25,115,64,13,0,128,64,21,0,-127,-191,-98,-177,-128,-255,-255,255,0,-128,-191,-223,64,-96,0,-127,191,32,46,105,-75,-64,-159,101,64,-96,-79,128,64,-96,71,-92,-174,213,-21,138,59,-255,255,0,-128,-191,-223,-239,-96,175,-32,143,32,-112,52,48,95,-80,-168,64,-96,175,-37,64,-96,-175,-58,-157,206,16,56,-16,-6,255,0,-128,-191,-223,-239,-247,175,40,96,32,-112,183,-31,104,-14,-135,195,-30,143,-40,114,-71,-163,-209,-151,200,17,119,36,110,-73,0,128,-64,-159,-207,231,12,-20,28,-16,102,40,36,30,-5,29,102,-12,-2,-21,-138,12,-14,134,-60,17,12,-106,-37,-29,-2,75,127,64,-96,-175,215,4,54,4,42,16,-72,52,-102,-69,73,-2,-22,133,9,-91,-173,103,39,-20,114,22,-46,151,50,34,-24,34,63,96,-80,167,4,-125,-58,-87,-87,-95,-163,-89,-172,-153,67,-78,125,61,-10,-116,186,4,20,-103,9,49,77,-34,5,-64,84,58,96,80,88,4,-125,190,-22,4,-63,-55,155,-43,149,9,50,-50,46,-97,-4,122,14,8,-44,-139,-108,-19,-59,142,17,-6,-43,24,80,-88,4,-125,190,32,46,46,62,47,20,72,18,-40,56,-66,-84,169,-1,-4,121,10,-21,138,-59,5,62,55,-8,61,36,-26,-88,171,-42,148,32,-111,52,-31,7,-20,-116,29,-113,-129,-50,-124,161,39,45,76,69,79,41,11,-94,84,20,-96,98,-23,7,102,171,42,107,32,-111,183,-31,9,-60,-48,151,-45,150,7,85,-57,32,-108,8,73,-91,0,-38,-120,-159,-8,-95,168,-2,-2,-82,16,42,-107,32,-111,183,36,-30,-48,-84,83,35,57,-34,-74,-13,94,-62,158,-17,-64,160,-41,0,128,64,53,56,39,25,55,-45,79,-107,181,-37,146,36,110,-40,-24,86,5,-94,-52,145,0,23,2,100,-37,5,105,40,48,46,-26,-15,22,-11,108,-5,21,-67,-33,181,37,109,36,110,-73,122,-25,37,-25,-140,133,-55,-72,77,-47,-35,145,5,48,-104,50,22,-117,-126,-82,119,68,-86,-77,129,54,37,-109,36,110,-73,-164,19,18,37,-32,-143,-36,145,27,25,51,75,-55,-6,-75,-165,-82,-33,144,-54,16,41,83,55,-23,-4,-2,-109,182,37,-73,-164,-209,-109,-68,89,-41,148,36,-50,-42,-69,63,-28,-141,26,-113,184,-30,36,-28,-109,-27,83,-74,40,47,36,-32,182,-37,-73,-164,-209,-232,180,23,-48,82,9,-68,152,-4,-88,-68,-122,-189,-75,162,18,-62,-64,-141,171,2,-41,148,-52,-39,-80,105,-37,146,-55,-155,-205,-230,-31,-48,83,23,96,30,46,23,-26,65,-95,-175,-5,18,118,-47,-5,-130,-8,-37,76,43,58,4,-42,43,146,55,-100,-178,-216,-236,-143,134,2,48,-66,61,97,-59,89,-34,-137,196,-23,92,-68,123,-45,-148,-127,116,5,105,13,75,74,-50,55,100,-77,-166,-211,233,11,-19,51,1,-118,63,76,27,20,35,100,-6,-5,-45,-150,14,-30,143,56,13,25,-75,-77,-25,-27,78,100,77,-89,-172,213,1,-40,-87,73,-54,-154,69,90,-40,-71,78,8,-123,45,-104,179,-31,47,7,-2,-20,-36,145,-55,95,23,1,77,89,-83,169,1,127,9,3,49,-36,-146,67,83,6,30,-43,-90,173,-7,104,23,50,-28,-105,102,21,55,45,-72,39,3,-55,89,83,86,1,127,-64,66,-37,6,-64,-146,96,63,28,24,-66,113,49,19,-53,154,39,32,-88,28,-28,-40,148,53,6,-67,-81,83,-86,1,127,-64,-160,5,58,16,-67,-161,59,-6,-17,12,30,49,-24,-47,127,16,-83,-58,-151,1,-50,103,-41,3,33,11,53,-86,170,42,-64,-160,-5,-46,-16,2,-33,-124,-4,-58,-29,-46,-9,-6,131,13,12,-113,20,27,8,51,-13,30,-10,45,-8,-7,43,170,-42,-64,-160,27,-104,-47,-45,-62,-113,7,-79,-10,-84,4,1,123,36,14,-98,2,94,8,98,-11,-30,-11,131,1,-61,29,76,-42,149,-53,27,-77,89,10,48,-59,11,122,0,55,-64,5,-54,19,27,19,14,-107,1,23,15,20,8,57,-2,32,-9,-33,-16,149,53,48,-104,76,90,0,14,61,54,3,39,1,18,6,9,17,3,6,-86,0,39,6,104,0,-11,0,62,0,-27,-8,38,53,-154,-104,76,90,-45,-17,45,9,-45,-44,0,-14,26,4,2,-7,16,7,0,-85,2,10,1,-14,1,28,-1,21,0,-11,-2];

const f0_he = [-255,-255,-255,-255,-255,-255,255,-16,135,-60,-157,-206,231,12,122,67,-94,175,-40,-148,201,27,-114,-185,220,18,-119,-187,221,-17,136,-59,-255,-255,-255,-255,-255,255,0,56,0,-114,-184,220,0,55,0,-14,134,0,-94,175,0,-44,-149,202,0,-51,-153,204,0,76,0,0,-255,-255,-255,-255,255,0,-48,12,64,-87,171,0,-61,-23,60,-20,-3,74,67,0,57,-8,132,0,59,-16,136,0,-68,6,0,2,-255,-255,-255,255,0,-64,-16,40,28,100,0,-61,17,23,-22,-11,61,-33,1,54,30,27,0,60,30,34,0,-68,-15,0,-19,-27,-255,-255,255,0,-64,-16,-32,23,-46,2,-59,15,-22,0,-5,-13,-24,-4,28,25,-42,-1,50,28,-45,1,-56,-11,0,-6,21,-6,-255,255,0,-127,0,-40,-20,-34,3,-37,11,-21,-1,-3,-10,12,-3,-9,11,-38,3,13,21,-44,-1,-16,-9,0,8,11,-2,10,255,0,-127,0,-32,-20,-26,4,-25,9,-19,0,-2,-11,3,-1,-6,7,-30,2,17,13,-40,-1,-20,-5,0,7,10,0,13,2,0,128,-64,-64,-16,30,5,-3,14,-15,0,8,12,1,-6,0,3,-26,-5,21,-15,-4,-7,23,-18,0,9,-10,-9,-1,-5,16,127,64,-95,-32,24,-3,5,-10,-2,-4,7,7,-3,-13,-1,4,-22,-7,14,-20,-4,8,20,21,0,0,-12,4,0,0,2,3,64,96,-80,48,-20,-11,-8,4,-17,10,-27,-5,-19,-6,2,-26,0,-17,-6,-5,11,-26,10,0,-5,12,11,-1,7,-3,11,8,96,80,88,-40,0,42,-6,-3,-24,-24,6,-20,1,-2,29,-1,-8,-60,-2,14,27,6,0,10,-29,1,-1,6,3,-1,-9,7,80,-88,4,0,43,8,4,16,-7,5,31,8,0,-21,5,-6,-52,0,4,36,8,0,5,38,-1,1,11,-9,-1,5,4,2,-88,171,-42,86,-10,-9,13,7,-9,29,-23,-2,-28,0,-11,-48,4,-21,39,9,0,-23,42,-10,-7,28,-1,-1,12,-4,11,12,171,42,107,-21,0,43,-8,-2,-29,-23,4,-29,-3,-2,42,0,-10,-65,3,0,25,-6,-5,6,-39,3,-3,11,5,0,-10,4,42,-107,32,0,45,9,-2,15,0,6,34,5,-1,-26,2,-5,-55,-2,0,31,8,7,0,50,3,3,16,-9,0,6,2,2,-107,181,-37,90,-9,9,-13,-8,8,23,-22,1,31,1,8,9,-1,0,15,8,-8,-2,41,10,10,4,-9,3,-14,7,0,-2,181,37,109,-18,0,-34,10,-4,-4,-18,-5,30,-3,1,12,-1,0,-13,3,-11,10,18,6,-13,-10,-9,3,-7,-5,-2,7,0,37,-109,36,0,-45,-14,6,9,-13,-3,-38,2,-7,-1,6,0,-11,3,-8,11,-28,-3,-17,14,-8,0,-14,6,9,-1,7,-1,-109,182,36,-91,0,-41,-9,-24,1,-31,9,-4,-1,1,0,4,4,-13,7,-28,5,14,12,-31,-4,-13,6,10,0,9,1,7,182,-36,-73,0,-36,-19,18,5,4,5,-20,3,4,0,1,-1,-12,-7,-12,5,-19,2,-24,5,10,2,2,6,-6,1,-5,5,-37,146,-55,-73,-14,-31,0,-18,7,-19,1,12,0,2,-2,6,-6,-8,-2,-21,-5,8,-4,-25,8,-8,10,13,-3,8,-5,9,146,55,-100,-27,-25,-24,-16,1,-20,-3,12,0,2,-7,2,-6,-5,9,-20,-7,14,17,-26,-10,-15,-3,16,8,8,2,9,1,55,100,-77,-50,-19,28,5,-2,17,11,0,3,9,2,-1,0,4,-24,-7,16,-19,-4,-13,19,22,1,9,-10,-11,-1,-6,15,100,77,-89,-39,22,-4,-3,-8,11,0,4,6,3,-13,-4,1,-20,5,13,-25,5,17,17,-26,-4,-7,-14,7,-1,1,5,-3,77,89,-83,44,-21,11,5,4,0,4,-22,4,19,2,2,5,4,17,-6,4,-19,-3,-14,6,23,2,4,0,-10,5,-1,1,89,83,86,-41,0,-32,5,-1,-1,4,-3,17,6,3,-4,14,6,14,-1,-19,-8,-2,-2,-28,7,-1,-1,3,8,10,-3,13,83,-86,2,0,-43,-13,-1,-2,3,-2,-33,-1,-1,-6,1,2,10,-1,-11,-4,2,-3,-24,-8,0,-2,-1,-3,0,1,0,-1,-86,171,42,-85,0,0,4,13,4,-12,-9,0,-4,-1,-1,4,-1,22,-3,3,-3,5,-9,0,3,-2,-3,-1,2,0,-1,-1,171,-42,-64,0,13,9,0,-4,-3,-5,0,-8,-2,-7,4,0,22,-1,1,-4,7,-1,0,-5,4,-1,-1,7,0,-1,1,1,-42,149,-53,26,12,0,-3,0,-5,0,6,3,9,0,0,-11,0,1,-2,14,0,0,-4,-6,-5,-1,0,2,-3,-2,-3,2,149,53,48,-53,13,29,-4,7,2,-3,0,4,-1,1,1,-1,1,4,1,4,-2,10,0,7,-1,-3,0,-1,1,1,-1,1,53,-154,-104,0,32,-8,-13,26,3,0,-7,-1,-2,-2,0,0,0,1,1,0,-7,0,3,1,8,0,-1,-1,-1,0,-3,3];

const f1_le = [-255,-255,-255,-255,-255,255,64,96,-80,-167,-211,233,11,122,66,-94,175,-40,-148,202,27,114,0,-128,18,-119,-187,-221,12,-122,-59,-94,-255,-255,-255,-255,255,0,16,0,-124,-189,222,0,56,0,-14,134,0,-94,175,0,44,0,-64,9,-55,-143,-199,8,-59,-44,-62,0,-255,-255,-255,255,0,-64,8,48,-102,178,0,-61,-23,60,-20,-4,74,67,0,-57,0,52,0,-35,-20,-126,3,-54,-4,-33,0,-6,-255,-255,255,0,-64,-16,16,20,113,0,-61,17,22,-22,-11,61,-33,1,-54,0,23,-2,-18,-13,-37,1,-49,-3,-12,0,-3,0,-255,255,0,-128,0,-40,6,-34,1,-60,17,-22,0,-5,-13,-24,-5,-27,0,-24,-1,14,-8,-24,0,-31,-4,-10,0,-2,0,-1,255,0,-128,0,-32,-20,-15,4,-40,10,-21,-1,-3,-10,12,-3,9,0,-20,0,6,5,-13,0,-9,-2,-9,0,0,0,-1,0,0,128,-64,-64,-16,30,-5,31,-23,-23,-2,10,-13,-2,13,13,0,17,-4,9,2,13,1,12,4,11,0,1,1,1,0,-2,128,64,-96,-32,24,-3,-21,-14,-30,-5,10,7,0,-27,5,1,-9,3,9,-2,-2,0,0,-3,-1,0,1,2,-2,0,0,0,64,96,-80,48,-20,-10,-4,-4,-27,12,-29,-6,-18,7,1,-14,0,-7,-4,2,0,0,-2,1,0,0,3,-2,0,0,0,0,96,80,88,-40,0,43,21,19,35,-32,7,-19,2,0,-20,0,-3,3,2,1,0,1,1,0,-1,4,-4,0,0,0,0,0,80,-88,4,0,43,8,-9,24,5,9,33,-8,-2,-18,0,-1,-5,0,0,0,0,0,0,1,0,-7,0,0,0,0,0,0,-88,171,-42,86,-10,-9,2,-7,-14,34,23,0,-7,-4,-1,-9,2,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,0,171,42,107,-21,0,43,35,10,43,26,0,-11,1,-4,-10,3,0,0,0,-1,2,0,-2,0,0,0,0,0,0,0,0,0,42,-107,32,0,45,9,-4,29,-15,0,14,-3,0,-9,0,0,-1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,-1,-107,181,-37,90,-9,9,-2,-15,0,30,-4,1,-2,3,0,4,-3,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,181,37,109,-19,0,-34,13,-1,21,-27,2,2,2,0,-2,-2,1,1,-1,0,0,0,1,0,0,0,0,0,0,0,0,0,37,-109,36,0,-45,14,0,7,9,-5,-13,0,-4,-6,2,1,8,1,-1,0,0,1,0,-2,0,0,0,0,0,0,0,0,-109,182,37,-91,0,-10,9,-14,5,-45,2,0,-1,1,-2,2,1,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,182,-37,-73,0,7,-32,-18,2,-7,9,-1,2,0,2,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-37,146,55,27,-22,28,7,-10,-10,-1,-4,0,1,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,146,-55,0,-39,22,2,-12,-5,1,-2,-1,1,0,1,0,0,0,1,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,-55,100,-77,-50,-42,-14,-4,-2,-8,0,-1,1,1,0,0,0,0,0,-1,0,0,0,0,0,0,0,1,0,1,0,0,0,100,-77,0,-42,-21,24,2,3,-1,-1,10,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,-77,0,-42,21,27,16,-1,0,1,1,-3,-1,1,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-42,21,11,16,-4,0,-3,-7,-3,-5,0,1,0,0,0,1,1,0,-1,0,0,-1,0,0,1,0,0,0,0,0,0,-83,0,-1,11,8,5,-11,20,-3,-15,0,1,-2,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,85,-11,24,6,-26,-1,-17,-5,-2,5,0,-1,0,0,0,1,0,-1,0,0,0,1,-1,0,0,0,0,0,0,0,0,170,42,0,27,-25,5,-40,-5,0,9,0,0,0,0,-1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,42,0,-27,-25,-3,41,-1,1,19,-10,4,0,0,1,0,0,-2,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,-27,12,-19,37,1,-1,6,1,5,3,-1,-2,-1,0,-1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,-53,-154,0,64,-3,9,0,-6,1,5,0,0,0,1,-1,0,0,0,0,0,0,1,0,0,0,0,2,0,0,0,0,0,-154,0,90,-3,14,20,9,-1,0,-4,0,-1,0,-1,-1,0,0,1,-1,0,0,0,-1,1,0,2,0,0,0,0,0,1];

const f1_he = [-255,-255,-255,-255,-255,255,64,96,-80,-167,-211,233,11,122,67,-94,175,-40,-148,201,27,114,0,-128,18,-119,-187,-221,12,-121,-59,-94,-255,-255,-255,-255,255,0,16,0,-124,-189,222,0,55,0,-14,134,0,-94,174,0,44,0,-64,9,-55,-143,-199,8,-59,-44,-62,0,-255,-255,-255,255,0,-64,8,48,-102,178,0,-61,-23,60,-20,-3,74,67,0,-57,0,52,0,-35,-20,-126,3,-54,-4,-33,0,-6,-255,-255,255,0,-64,-16,16,20,113,0,-61,17,22,-22,-11,61,-33,1,-54,0,23,-2,-18,-13,-37,1,-49,-3,-12,0,-3,0,-255,255,0,-127,0,-40,6,-34,2,-60,17,-22,0,-5,-13,-24,-4,-27,0,-24,-1,14,-8,-24,0,-31,-4,-10,0,-2,0,-2,255,0,-127,0,-32,-20,-15,5,-40,10,-21,-1,-3,-11,12,-2,10,0,-20,0,6,5,-13,0,-9,-3,-9,0,-1,0,-1,0,0,127,-64,-64,-16,30,-5,31,-23,-23,-2,11,-13,-2,13,13,0,17,-4,9,2,13,1,12,4,11,0,1,1,1,0,-2,128,64,-96,-32,24,-3,-21,-14,-30,-5,10,7,0,-27,5,0,-14,2,9,-4,-3,0,0,-4,1,0,1,2,-2,0,0,0,64,96,-80,48,-20,-10,-4,-4,-27,12,-29,-5,-18,7,0,-16,0,-6,-6,0,0,0,-2,0,0,0,3,-2,0,0,0,0,96,80,88,-40,0,43,21,19,35,-32,7,-19,1,0,-16,-1,-5,5,0,0,0,2,0,0,0,3,-4,0,0,0,0,0,80,-88,4,0,43,8,-9,23,5,9,33,-8,0,-18,0,-3,-4,0,0,0,0,0,2,1,0,-8,0,0,0,0,0,0,-88,171,-42,86,-10,-9,2,-7,-14,34,24,0,-8,0,-2,-5,1,0,-1,0,0,-1,0,0,0,0,0,0,0,0,0,0,171,42,106,-21,0,43,35,10,43,26,0,-10,0,-4,-5,0,0,1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,42,-107,32,0,45,9,-3,28,-15,0,14,1,0,-6,1,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,-107,181,-37,90,-9,9,-2,-14,0,18,-1,0,-3,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,181,37,109,-19,0,-34,13,0,11,-1,-1,4,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,37,-109,36,0,-45,14,0,8,6,-1,-6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-109,182,37,-91,0,0,17,-10,3,-8,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,182,-37,-73,0,0,-13,-6,3,-5,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-36,146,55,0,-7,-9,5,-12,-5,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,146,-55,-1,-14,-6,-15,-8,-5,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-55,0,-39,-10,-14,13,-2,1,-2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,-39,-22,-15,11,6,1,2,3,0,2,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,-39,-22,-31,14,6,1,0,0,-2,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-22,-31,27,1,-2,-3,-3,0,-2,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-32,26,-1,-7,1,-7,-4,0,-3,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,28,57,0,4,-3,-7,1,-5,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,57,28,-11,-4,-11,-4,-5,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,28,-24,-22,-12,-6,22,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-24,-22,6,-22,22,-2,-1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-22,-110,-15,57,-1,4,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-109,30,77,-3,3,-4,2,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1];

const f2_le = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1];

const f2_he = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1];

const f3_le = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

const f3_he = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

viz_bias('g0_le', g0_le, 32,32);
viz_bias('g0_he', g0_he, 32,32);
viz_bias('f0_le', f0_le, 32,32);
viz_bias('f0_he', f0_he, 32,32);
viz_bias('f1_le', f1_le, 32,32);
viz_bias('f1_he', f1_he, 32,32);
viz_bias('f2_le', f2_le, 32,32);
viz_bias('f2_he', f2_he, 32,32);
viz_bias('f3_le', f3_le, 32,32);
viz_bias('f3_he', f3_he, 32,32);

//var leg_32 = viz_legend('heat_map',2,32); 
</script>
