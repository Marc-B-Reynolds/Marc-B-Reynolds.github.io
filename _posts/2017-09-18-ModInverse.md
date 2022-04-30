---
layout:       post
title:        Integer multiplicative inverse via Newton's method
categories:   [math]
tags:         [integer]
description:  A short note on computing the modular multiplicative inverse of an odd integer.
---

\\
I wanted to mention computing the multiplicative inverse via Newton's iteration in a post and to my suprise I couldn't find a simple overview which is publically available.  (We interrupt this post in final edit with just such a thing from Daniel Lemire[^lemire]).  For a formal presentation this 2012 paper[^newton]  by Jean-Guillaume Dumas is pretty concise and easy to understand. Note I'll use mod as an operator instead of for congruence and will only talk about power-of-two modulus (linked paper covers generalization). First let's define the function.  Given an odd integer $a$ we want to find integer $x$ such that:

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
static inline uint32_t mod_inverse_3(uint32_t a)
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
Notice if we instead had 4 bits for our initial value we'd have after each successive step: $\left(8,16,32,64,128,\ldots\right)$ bits so we'd need only three iterations.  So if we can get an extra bit for cheaper than the cost of one step then we're ahead. There are two quadradic solutions[^quadradic] :

$$ \begin{eqnarray}
a^{-1} \bmod 16 & = & \left(a + a^2 - 1\right) \bmod 16 \\
                & = & \left(a - a^2 + 1\right) \bmod 16
\end{eqnarray}$$

\\
Using the first one we now have:

{% highlight c %}
static inline uint32_t mod_inverse_4(uint32_t a)
{
  uint32_t x;
  x = (a*a)+a-1;     //  4 bits (For serial comment below: a*a & a-1 are independent) 
  x *= 2-a*x;        //  8
  x *= 2-a*x;        // 16
  x *= 2-a*x;        // 32
  return x;
}
{% endhighlight %}

\\
Another paper[^mmul] mentions a method for initial value good to 5 bits:

{% highlight c %}
static inline uint32_t mod_inverse_5(uint32_t a)
{
  uint32_t x;
  x = 3*a ^ 2;       //  5 bits
  x *= 2-a*x;        // 10
  x *= 2-a*x;        // 20
  x *= 2-a*x;        // 40 -- 32 low bits
  return x;
}
{% endhighlight %}


\\
Dumas carries through with the derivation (SEE: Section 3.3) to produce algorithm 3, which reworked looks like:

{% highlight c %}
static inline uint32_t mod_inverse_d(uint32_t a)
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
The first three variants are almost identical in structure so it looks[^looks] like an architectural toss up on performance.  The last I find quite interesting.  The other two require every operation to be performed in serial order...not so for the third.  Renamed and with ops that can happen at the same time on the same line.

     1: u  = 2-a    i0 = a-1
     2: i1 = i0*i0
     3: i2 = i1*i1  t1 = i1+1
     4: i3 = i2*i2  t2 = i2+1  u *= t1
     5: i4 = i3*i3  t3 = i3+1  u *= t2
     6: t4 = i4+1              u *= t3
     7: u *= t4
     (done)
     
     
\\
The Dumas version cannot drop any of the iteration steps, but previous methods can.  This is potentially interesting if there's a known bit-bound on input values.  We have choices between 3, 4 and 5 inital bits and must perform at least one step.

<br>

------

#### Jeffrey Hurchalla's method <small>ADDED: 20220407</small>

\\
Let's take stock of where we are. The classic version has the downside of being a serial computation but we we can choose the method of inital value which in turn controls how many iteration we need to perform. As far as I know the best (in terms of cost) "computed" intial value is the 5-bit of Montgomery. For Dumas we lose the initial value feature and in trade get some instruction level parallelism.

In April of 2022 Jeffrey Hurchalla released a paper[^hurchalla] which gives us both: choose inital value and broken dependency chains (see paper for details). Here's an example implementation for 32-bit integers using the 5-bit initial value:

{% highlight c %}
static inline uint32_t mod_inverse_h(uint32_t a)
{
  uint32_t x = (3*a)^2; 
  uint32_t y  = 1 - a*x;
  x = x*(1+y); y *= y;
  x = x*(1+y); y *= y;
  x = x*(1+y);
  return x3;
}
{% endhighlight %}


\\
The original version of this post didn't talk about performance at all. All of these are going to be very close for a general purpose implementation and the "best" is going to depend on the exact micro-architecture and surrounding code. As a rule of thumb the Hurchualla variant is the most promising.


<br>

------

#### Spitballing performance <small>ADDED: 20220407</small>

\\
All of these functions are small and fast in an absolute sense. So the only reason to look closely at various alternate choices is the case of computing tons of them when there is very little other work happening at the same time. So it's probably a waste of time to read any further.

But if you're still here let's look at why I'm suggest Hurchualla as the general purpose choice. I'm using the intel family of processors for following:

<div class="container">
  <div class="row">
    <div class="col-sm">

{% highlight nasm %}

; mod_inverse_5 
; (5-bit, serial)
lea    ecx, [rdi+2*rdi]
xor    ecx, 2
mov    edx, ecx
imul   edx, edi
mov    eax, 2
mov    esi, 2
sub    esi, edx
imul   esi, ecx
mov    ecx, esi
imul   ecx, edi
mov    edx, 2
sub    edx, ecx
imul   edx, esi
imul   edi, edx
sub    eax, edi
imul   eax, edx
{% endhighlight %}

</div>
<div class="col-sm" markdown="1">

{% highlight nasm %}

; mod_inverse_h 
; (5-bit, Hurchalla)
lea    ecx, [rdi+2*rdi]
xor    ecx, 2
imul   edi, ecx
mov    eax, 1
sub    eax, edi
mov    edx, 2
sub    edx, edi
imul   edx, ecx
imul   eax, eax
lea    ecx, [rax+1]
imul   ecx, edx
imul   eax, eax
inc    eax
imul   eax, ecx
{% endhighlight %}

</div>
<div class="col-sm" markdown="1">

{% highlight nasm %}

; mod_inverse_d
; (Dumas)
mov    eax, 2
sub    eax, edi
dec    edi
imul   edi, edi
lea    ecx, [rdi+1]
imul   ecx, eax
imul   edi, edi
lea    edx, [rdi+1]
imul   edx, ecx
imul   edi, edi
lea    eax, [rdi+1]
imul   eax, edx
imul   edi, edi
inc    edi
imul   eax, edi
{% endhighlight %}

</div>
</div>
</div>

\\
In terms of both high level opcodes they are all close:


{: .center }
|               | mod_inverse_5 | mod_inverse_h | mod_inverse_d | notes |
| :---:         | :---:         | :---:         | :---:         | ---  |
| total ops     | 17            | 14            | 15            | mirco-op count is the same    |
| LEA           |  1            |  2            |  3            | latency of 1, ports {1,5}     |
| IMUL          |  6            |  6            |  8            | latency of 3, ports {1}       |
| other         |  9            |  6            |  4            | latency of 1, ports {0,1,5,6} |

\\
OK. I really should have use 64-bit (each need one extra round) here but I'm too lazy to go back in modify. But anyway Hurchalla has the fewest number of uOps to issue. Compared to Dumas (only one more uOP) we're eliminating 2 multiplies (latency 3, tied to port 1) and a LEA. Few of Dumas ops can go through ports {0,1,5,6} adding to the probability of competition with surrounding code. Compared with the "classic" serial method: the classic needs to execute 3 more uOps, has one less LEA then Hurchalla and 3 more uOps can go through any of {0,1,5,6} (the same number of extra uOps to be issued). Bottom line is that any of these could end up being "fastest" at some specific site but Hurchalla on average will win out.

Now back to the hypothetical problem of computing many inverses with little other work happening at the same time. In that case we can very likely afford to use a small look-up table. With a 256 byte table we can drop one round giving:

<div class="container">
  <div class="row">
    <div class="col-sm">

{% highlight c %}
// Hurchalla with 8-bit initial value from table
inline uint32_t mod_inverse_th(uint32_t a)
{
  uint32_t x = (uint32_t)(mod_inverse_table_8[a & 0xff]);
  uint32_t y = 1-a*x;
  x = x*(1+y); y *= y;
  x = x*(1+y);
  return x;
}
{% endhighlight %}

<div class="card"><div class="card-header"><details markdown="1">
<summary>the data table</summary>

{% highlight c %}
// 8-bit mod inverse table. only odd
// elements accessed for legal input
const uint8_t mod_inverse_table_8[] =
{
  0,0x01,0,0xab,0,0xcd,0,0xb7,
  0,0x39,0,0xa3,0,0xc5,0,0xef,
  0,0xf1,0,0x1b,0,0x3d,0,0xa7,
  0,0x29,0,0x13,0,0x35,0,0xdf,
  0,0xe1,0,0x8b,0,0xad,0,0x97,
  0,0x19,0,0x83,0,0xa5,0,0xcf,
  0,0xd1,0,0xfb,0,0x1d,0,0x87,
  0,0x09,0,0xf3,0,0x15,0,0xbf,
  0,0xc1,0,0x6b,0,0x8d,0,0x77,
  0,0xf9,0,0x63,0,0x85,0,0xaf,
  0,0xb1,0,0xdb,0,0xfd,0,0x67,
  0,0xe9,0,0xd3,0,0xf5,0,0x9f,
  0,0xa1,0,0x4b,0,0x6d,0,0x57,
  0,0xd9,0,0x43,0,0x65,0,0x8f,
  0,0x91,0,0xbb,0,0xdd,0,0x47,
  0,0xc9,0,0xb3,0,0xd5,0,0x7f,
  0,0x81,0,0x2b,0,0x4d,0,0x37,
  0,0xb9,0,0x23,0,0x45,0,0x6f,
  0,0x71,0,0x9b,0,0xbd,0,0x27,
  0,0xa9,0,0x93,0,0xb5,0,0x5f,
  0,0x61,0,0x0b,0,0x2d,0,0x17,
  0,0x99,0,0x03,0,0x25,0,0x4f,
  0,0x51,0,0x7b,0,0x9d,0,0x07,
  0,0x89,0,0x73,0,0x95,0,0x3f,
  0,0x41,0,0xeb,0,0x0d,0,0xf7,
  0,0x79,0,0xe3,0,0x05,0,0x2f,
  0,0x31,0,0x5b,0,0x7d,0,0xe7,
  0,0x69,0,0x53,0,0x75,0,0x1f,
  0,0x21,0,0xcb,0,0xed,0,0xd7,
  0,0x59,0,0xc3,0,0xe5,0,0x0f,
  0,0x11,0,0x3b,0,0x5d,0,0xc7,
  0,0x49,0,0x33,0,0x55,0,0xff
};
{% endhighlight %}

</details></div></div>
</div>
<div class="col-sm" markdown="1">

{% highlight nasm %}
movzx  eax, dil
movzx  ecx, byte ptr [rax+mod_inverse_table_8]
imul   edi, ecx
mov    eax, 1
sub    eax, edi
mov    edx, 2
sub    edx, edi
imul   edx, ecx
imul   eax, eax
inc    eax
imul   eax, edx
{% endhighlight %}

</div>
</div>
</div>

<br>

<div class="container">
  <div class="row">
    <div class="col-sm">

{% highlight c %}
// classic with 8-bit initial value from table
uint32_t mod_inverse_t5(uint32_t a)
{
  uint32_t x = (uint32_t)(mod_inverse_table_8[a & 0xff]);
  x *= 2-a*x;
  x *= 2-a*x;
  return x;
}
{% endhighlight %}

</div>
<div class="col-sm" markdown="1">

{% highlight nasm %}
movzx  eax, dil
movzx  ecx, byte ptr [rax+mod_inverse_table_8]
mov    edx, ecx
imul   edx, edi
mov    eax, 2
mov    esi, 2
sub    esi, edx
imul   esi, ecx
imul   edi, esi
sub    eax, edi
imul   eax, esi
{% endhighlight %}
</div>
</div>
</div>

<br>

\\
No clear winner here between the two. We drop one round for a data dependent load (5 cycle latency & issued to port 2 or 3). Either should outperform the computed initial value variants.


------


References and Footnotes
------

[^nrmethod]:   **"Wikipedia: Newton's method"**, ([page](http://en.wikipedia.org/wiki/Newton%27s_method))
[^newton]:     **"On Newton-Raphson iteration for multiplicative inverses modulo prime powers"**, Jean-Guillaume Dumas 2012 ([PDF](http://arxiv.org/pdf/1209.6626v2.pdf))
[^hensel]:     **"Wikipdia: Hensel's lemma"**   ([page](https://en.wikipedia.org/wiki/Hensel%27s_lemma))
[^lemire]:     **"Computing the inverse of odd integers"**, Daniel Lemire 2017 ([page](http://lemire.me/blog/2017/09/18/computing-the-inverse-of-odd-integers/))
[^looks]:      AKA I've spent zero time really thinking about performance here.
[^quadradic]:  I just brute forced this in Mathematica.
[^mmul]:       **"Efficient long division via Montgomery multiply"**, Ernst W. Mayer, 2016 ([PDF](http://arxiv.org/abs/1303.0328))
[^hurchalla]:  **"Speeding up the Integer Multiplicative Inverse (modulo 2<sup>w</sup>) "**, Jeffrey Hurchalla, 2022 ([PDF](https://arxiv.org/abs/2204.04342))












