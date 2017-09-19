---
layout:       post
title:        Integer multiplicative inverse via Newton's method
categories:   [math]
tags:         [integer]
description:  A short note on computing the modular multiplicative inverse of an odd integer.
---

\\
I wanted to mention computing the multiplicative inverse via Newton's iteration in a post and to my suprise I couldn't find a simple overview which is publically available.  (We interrupt this post in final edit with just such a thing from Daniel Lemire[^lemire]).  For a formal presentation this 2016 paper[^newton]  by Jean-Guillaume Dumas is pretty concise and easy to understand. Note I'll use mod as an operator instead of for congruence and will only talk about power-of-two modulus (linked paper covers generalization). First let's define the function.  Given an odd integer $a$ we want to find integer $x$ such that:

$$ \begin{equation}
\left(a~x\right) \bmod 2^\beta = 1
\end{equation} $$

\\
For reals and floating point we can compute the inverse of $a$ given a good inital guess $x_0$ and applying the following recurrence (Newton's method[^nrmethod]) a sufficient number of times:

$$ \begin{equation}
x_{n+1} = x_n\left(2-a~x_n\right) \label{rr}
\end{equation} $$

\\
Very roughly speaking *Hensel's (lifting) lemma*[^hensel] extends Newton's method to modular arithmetic where (for our case) if we have the value $x_n$ such that:

$$ \left(a~x_n\right) \bmod 2^k = 1 $$

\\
then after applying the recurrence $\eqref{rr}$ we'd have:

$$ \left(a~x_{n+1}\right) \bmod 2^{2k} = 1 $$

\\
Which each application the number of bits double.  We need to compute an initial value $x_0$ that's correct to some $k$ bits. It happens that $a$ is its own inverse up to three bits $\left(k=3\right)$ so one option is to simply set $x_0 = a$.  Starting with 3-bits and repeated application of the recurrence would give: $\left(6,12,24,48,96,\ldots\right)$.  So assuming working in 32-bit register with any odd 32-bit input we'd need four iterations:


{% highlight c %}
static inline uint32_t mod_inverse_1(uint32_t a)
{
  uint32_t x,t;
  x = a;             //  3 bits
  x *= 2-a*x;        //  6
  x *= 2-a*x;        // 12
  x *= 2-a*x;        // 24
  x *= 2-a*x;        // 48 / retaining bottom 32
  return x;
}
{% endhighlight %}

\\
Notice if we instead had 4 bits for our initial value we'd have after each successive step: $\left(8,16,32,64,128,\ldots\right)$ bits so we'd need only three iterations.  So if we can get an extra bit for cheaper than the cost of one step then we're ahead. There are two quadradic solutions:

$$
a + a^2 - 1 \\
a - a^2 + 1
$$

\\
Using the first one we now have:

{% highlight c %}
static inline uint32_t mod_inverse_2(uint32_t a)
{
  uint32_t x,t;
  x = (a*a)+a-1;     //  4 bits (For serial comment below: a*a & a-1 are independent) 
  x *= 2-a*x;        //  8
  x *= 2-a*x;        // 16
  x *= 2-a*x;        // 32
  return x;
}
{% endhighlight %}


\\
Dumas carries through with the derivation (SEE: Section 3.3) to produce algorithm 3, which reworked looks like:

{% highlight c %}
static inline uint32_t mod_inverse_3(uint32_t a)
{
  uint32_t u = 2-a;
  uint32_t i = a-1;
  i *= i; u *= i+1;
  i *= i; u *= i+1;
  i *= i; u *= i+1;
  i *= i; u *= i+1;
  return u;
}
{% endhighlight %}

\\
The first two variants are almost identical in structure so it looks[^looks] like an architectural toss up on performance.  The third I find quite interesting.  The other two require every operation to be performed in serial order...not so for the third.  Renamed and with ops that can happen at the same time on the same line.

     1: u  = 2-a    i0 = a-1
     2: i1 = i0*i0
     3: i2 = i1*i1  t1 = i1+1
     4: i3 = i2*i2  t2 = i2+1  u *= t1
     5: i4 = i3*i3  t3 = i3+1  u *= t2
     6: t4 = i4+1              u *= t3
     7: u *= t4
     (done)



<br>

------


References and Footnotes
------

[^nrmethod]:   *"Wikipedia: Newton's method"*, ([page](http://en.wikipedia.org/wiki/Newton%27s_method))
[^newton]:     *"On Newton-Raphson iteration for multiplicative inverses modulo prime powers"*, Jean-Guillaume Dumas 2016 ([PDF](http://arxiv.org/pdf/1209.6626v2.pdf))
[^hensel]:     *"Wikipdia: Hensel's lemma"*   ([page](https://en.wikipedia.org/wiki/Hensel%27s_lemma))
[^lemire]:     *"Computing the inverse of odd integers"*, Daniel Lemire 2017 ([page](http://lemire.me/blog/2017/09/18/computing-the-inverse-of-odd-integers/))
[^looks]:      AKA I've spent zero time really thinking about performance here.











