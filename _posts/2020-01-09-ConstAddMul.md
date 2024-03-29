---
layout:       post
title:        'FMA: extended precision addition and multiplication by constant'
categories:   [math]
tags:         [FMA,floating point]
description:  'Brushstrokes using FMAs to compute: K+x and Kx with K in higher than native precision.'
plotly:       false
---

\\
This is a quick note on a pair of small things that you either already know or will probably never need. OK, enough of my hard sell. As usually I'll use single precision as the working example and I'll only use *binary32* or *float* to talk about the FP format so I can use single to mean **one** of something.

We have a floating point value $x$ to which we want to either add or multiply by some constant $K$ to get a floating point value in the same format as $x$. But $K$ doesn't fit in that format and we want to bump up the accuracy of the result. As an example: binary32 $x$ and multiplying by $\pi$ the way normal people do:

{% highlight c %}
inline float f32_mul_pi(float x) { return ((float)M_PI)*x; }
{% endhighlight %}

\\
but we're nerding out and instead want a drop in replacement that is a little more accurate. Note this example is taking a double constant `M_PI` and casting to a *float* (which rounds it to a binary32). This works out in this case but in general a constant rounded to a double and then to a float is not the same as the constant directly rounded to a float.

We're going to need some extra tools to work with and I'll write using the following freely available: 

* [Sollya](http://sollya.gforge.inria.fr/) is a free library and DSL designed for building math functions.
* The [WolframAlpha](https://www.wolframalpha.com/) website to perform prime factorization of integers.

\\
Recall that *fused multiply add* (FMA) operations (exposed in C/C++ as `fma` and `fmaf`) computes `a*b+c` with a single rounding step:

$$ \text{RN}\left(a \cdot b + c\right) $$


\\
where $\text{RN} \left(\cdot\right)$ is simply the hardware *round to nearest (ties to even)* for the given format.

<br>

<div class="alert alert-info" role="alert" markdown="1">
**WARNING:** Only recent updates of Visual C (2019) treat `fma` and `fmaf` as intrinsics and older versions will instead calls some obnoxiously expensive functions.
</div>

<br>

------

Multiply by a constant <small>1 FMA + 1 product</small>
------

\\
 We can represent a number $\left(v\right)$ using more than one floating point value. I will call using the logical sum of two $\left(v = h+l\right)$ an *unevaluated pair*. We can build all kinds of tools out of them, but we need just one: multiply a pair with a single value and return a single value. A previous (and related) blog post ([*"Floating point division with known divisor"*](http://marc-b-reynolds.github.io/math/2019/03/12/FpDiv.html)) uses this and has a few more details.
 
{% highlight c %}
typedef struct { float h,l; } f32_pair_t;

// multiply 'x' by the value represented by the unevaluated pair p=(h,l)
inline float f32_up_mul(const f32_pair_t* const p, float x) 
{
  return fmaf(x, p->h, x*p->l);
}

{% endhighlight %}
  

\\
So the multiply by constant of this post is just a special case of this. We need to factor our constant $\left(K\right)$ into a pair:

 $$ \begin{align*}
 H & = \text{RN}\left(K\right) \\
 L & = \text{RN}\left(K-H\right)
 \end{align*} $$

\\
So $H$ is just the constant rounded to a single, $L$ is the error rounded to a single and $K \approx H+L$. Then we compute as:


$$ Kx \approx \text{RN}\left(H \cdot x + \text{RN}\left(L \cdot x\right)\right) $$

\\
Let's fire up sollya (`rlwrap sollya --nocolor`) and use it as a REPL to factor $\pi$:

    > single(Pi);
	3.1415927410125732421875
    > single(Pi-single(Pi));
	-8.74227765734758577309548854827880859375e-8
    > quit;

\\
(the function `single` rounds to *binary32* to nearest) and shoving the output into the above defs:

{% highlight c %}
static const f32_pair_t pi = {.h = 3.1415927410125732421875f, .l=-8.74227765734758577309548854827880859375e-8f};

inline float f32_mul_pi(float x) { return f32_up_mul(&pi, x); }
{% endhighlight %}


\\
So we have our drop-in replacement for the function from the introduction. What's interesting is that our replacement returns:

$$
\text{RN}\left(\pi~x\right)
$$

\\
which in English is: multiplies $x$ by $\pi$ as if in infinite precision and then rounds to nearest. The original version returns this ideal result about 66.8053% of the time. For a given constant and floating point format precision the FMA implementation will either: always return the correctly rounded result (like this $\pi$ example), always correctly rounded except for one input bit pattern (per power-of-two interval like the "division by constant"), or will be incorrect for multiple bit patterns. The methods to determine what happens for a given constant and FP format are described in [*"Correctly rounded multiplication by arbitrary precision constants"*](http://perso.ens-lyon.fr/jean-michel.muller/MultConstant.html)

Here's a sollya def to split out these constants as per this post:

    // computes a non overlapping pair out of input v
    procedure f32_mul_k(name,v)
    {
      var h,l,e;
      h = single(v);
      l = single(v-h);
      print("const f32_pair_t f32_mul_k_" @ name @
        " = {.h = " @ h @ "f, .l=" @ l @ "f};");
    };



\\
For your amusement here's a small table of some constants $\left(K \approx H+L\right)$, $CR$ indicates if always correctly round, ? = I don't know (I haven't implemented the testing methods) and '%' is percent of not correctly rounded when computed the normal way.

{: .center }
|K           | H             | L              | CR | % |
|:---:       | ---:          |---:            |:---: | ---:|
|$\pi$       |0x1.921fb6p1f  |-0x1.777a5cp-24f|Y | 33.194710 |
|$1/\pi$     |0x1.45f306p-2f | 0x1.b93910p-27f|Y | 48.123135 |
|$\log(2)$   |0x1.62e430p-1f |-0x1.05c610p-29f|Y |  3.260410 |
|$1/\log(2)$ |0x1.715476p0f  | 0x1.4ae0c0p-26f|Y | 15.840387 |
|$\log(10)$  |0x1.26bb1cp1f  |-0x1.12aabap-25f|Y | 16.824018 |
|$1/\log(10)$|0x1.bcb7b2p-2f |-0x1.5b235ep-27f|Y | 28.183519 |
|$e$         |0x1.5bf0a8p1f  | 0x1.628aeep-24f|? | 36.054657?|
|$1/e$       |0x1.78b564p-2f |-0x1.3a621ap-27f|? | 29.529118?|

\\
Something I haven't gotten around to explicitly mentioning is that since scaling by a power-of-two is error free that these precomputed unevaluated pairs can be simply scaled as well.

<br>

------

Carry it through...  <small>added: 20221022</small>
------

<br>

<div class="alert alert-info" role="alert" markdown="1">
**WARNING:** All of these are just examples. A good factorization depends on the range of input(s) and relative magnitudes of the individual quantities.
</div>


The previous by itself is mostly useless and just tightens up the final bit. However the real intent was for carrying through to rest of the computation.  Let's look at adding to it the subresult $c$ (the interesting case is when $c$ is significantly smaller than $Kx$):

$$ \begin{align*}
Kx + t &= \left(H+L\right)x+c \\
       &= Hx + Lx + c         \\
	   &= \left(Hx + \left(Lx + c\right)\right)
\end{align*} $$

\\
which is simply a nested FMA structure:

{% highlight c %}
// compute a*b+c where 'a' is represented by an unevaluated pair
static inline float f32_up_madd(const f32_pair_t* const a, float b, float c)
{
  return fmaf(a->h, b, fmaf(a->l, b, c));
}
{% endhighlight %}


\\
A more practical example usage when $c$ is a function of $x$ and we're computing a polynomial which is a minimax approximating of a function on some interval. Let's look at a polynomial in Horner's form (and without a constant term) example:

$$ \begin{array}{c}
x(c_1 + x(c_2 + \cdots ))  \\
x((H+L) + x(c_2 + \cdots )) \\
H~x + L~x + x^2(c_2 + \cdots ) \\
H~x + x(L + x~\underbrace{(c_2 + \cdots ))}_r)
\end{array} $$

\\
where $c_1$ is the unevaluated pair $H+L$, distribute the outer $x$, nest $L$ and denote a *tail* term $r$. So we have:

$$
H~x + x(r~x + L)
$$

which in translated into code give:

{% highlight c %}
  r = ...                      // previous terms like: 'r = fma(r,x,c2)' 
//r = x*fmaf(r,x,c1)           // straight version: vs.
  r = fmaf(H,x, x*fma(r,x,L))  //   pair version: +1 constant load, +1 fma
{% endhighlight %}

\\
so this expansion has the same cost as adding one more term to our approximation.

If our polynomial is *even* then we've precomputed $s=x^2$ and could just compute as the previous but let's consider another option. Replacing $s$ with $x^2$ and and reworking:

$$ \begin{array}{c}
x^2(c_1 + x^2(c_2 + \cdots ))  \\
s(c_1 + s(c_2 + \cdots ))  \\
\cdots \\
H~s + s(r~s + L) \\
H~x^2 + x^2(r~x^2 + L) \\
x(H~x + x(r~x^2 + L)) \\
x(H~x + x(r~s + L)) 
\end{array} $$

and in code:

{% highlight c %}
  r = ...                         // previous terms
//r = s*fmaf(r,s,c1);             // straight version: vs.
  r = fmaf(H,s, s*fma(r,s,L));    //   pair version 's':  +1 constant load, +1 fma  (as previous expansion)
   or
  r = x*fmaf(H,x, x*fmaf(r,s,L)); //   pair version 'x':  +1 constant load, +1 fma, +1 mul
{% endhighlight %}

\\
so choosing this second form adds more more product. And let's rinse and repeat for *odd* polynomials:

$$ \begin{array}{c}
x(c_1 + x^2(c_2 + \cdots ))  \\
x(c_1 + s(c_2 + \cdots ))  \\
\cdots \\
H~x + L~x + x^3 r  \\
H~x + x(L + x^2 r) \\
H~x + x(s~r + L) 
\end{array} $$

and again in code:

{% highlight c %}
  r = ...                        // previous terms & we've precomputed t=s*x
//r = fmaf(c1,x,t*r);            // straight version (1): x*c1 + x^3*r     (t=x^3 & computed earlier)
//r = x*fmaf(s,r,c1);            // straight version (2): x*(x^2*r + c1)
  r = fmaf(H,x, x*fmaf(r,s,L));  //   pair version: vs. (1) +1 constant load, +1 fma, -1 mul (t)
                                 //               : vs. (2) +1 constant load, +1 fma
{% endhighlight %}

which works out to be either an additional constant load (count fma & mul as equivalent in cost) or the cost of an extra term.

<br>

------

Overlapping pairs <small>added: 20221022</small>
------

\\
The constant factorization method given above produces a *non overlapping* pair which maximizes the precision of the constant. It can be desirable to produce a less accurate factorization. Here's a solla script that produces a pair where they both have the same sign.

    // produces same signed constants. generally
    // loses one bit of accuracy
    procedure f32_mul_samesign_k(name,v)
    {
      var ha,la;
      ha = round(v,24,RZ);
      la = single(v-ha);
    
      print("const f32_pair_t f32_mul_k_" @ name @
        " = {.h = " @ ha @ "f, .l=" @ la @ "f};");
    };

\\
The motivation here is when wanting to handle returning a signed zero for "free".  Given: `fma(a,b,c)` where `a*b` is zero and `c` is zero then the result is `-0` if both are `-0` and otherwise postive zero.

Another example is producing a $H$ with less than the maximum number of bit. FWIW: I've never used this...but I've never really tried either. But anyway this shrinks the magnitude gap between $H$ and $L$ (roughtly doubles for each bit) and the product $H*x$ requires fewer bits to be represented exactly. Why I think this could be interesting is because if we're bothering with extended precision for a polynomial approximaion then our coefficients are more than sufficient to hit a correctly rounded result and the problem has become managing the errors in performing the computation. The downside is we're losing some effectiveness of the *fma* with the $H$ term. 

    // produces constants where H is a (24-g) bits
    procedure f32_mul_p_k(name,v,g)
    {
      var ha,la;
      ha = round(v,24-g,RN);
      la = single(v-ha);
  
      print("const f32_pair_t f32_mul_k_" @ name @
        " = {.h = " @ ha @ "f, .l=" @ la @ "f};");
    };


<br>

------

Add a constant <small>1 FMA</small>
------

\\
Let me first note that this method **does not** help with [catastrophic cancellation](https://en.wikipedia.org/wiki/Loss_of_significance) since the source of that is the rounding error already present in $x$.

Like with the extended precision product this can be used as part of function approximations. Here's a [motivational example](https://forums.developer.nvidia.com/t/more-accurate-and-somewhat-faster-implementation-of-atanf/214453) implemenation of *atan* written by Norbert Juffa which needs to add $\frac{\pi}{2}$.

We can add by an extended precision constant as follows:

$$ K + x \approx \text{RN}\left(A \cdot B + x\right) $$

\\
So we need to again factor $K$ but this time into a multiplicative pair $\left(K \approx A \cdot B \right)$ instead of a sum. Here's some sollya defs to get started:

{% highlight c %}
procedure add_k(k,bits)
{
  var f,d,p,e,m,t;
  f = round(k,bits,RN);
  d = mantissa(round(k,bits,RU)-f);
  m = mantissa(f);
  e = exponent(f);
  p = precision(m);
  m = m*2^(bits-p);
  e = e-bits+p;

  return [|f,m,e,d|];
};

procedure add_k_q(k,bits) { var f; f=round(k,bits,RN); return [|mantissa(f),exponent(f)|]; };

procedure f32_add_k(x) { add_k(x, 48); };
procedure f64_add_k(x) { add_k(x,106); };

f32_add_k(Pi);
f64_add_k(Pi);
f32_add_k(2/(sqrt(5)+1));
add_k_q(2/(sqrt(5)+1),48);
{% endhighlight %}

\\
Breaking down `add_k`:

* `f` the constant rounded to nearest at twice the precision of the format (48/106 for float/double)
* `d` denotes if `f` was rounded up (0) or rounded down (1). This could be needed to break a 'tie'.
* `m` significand (or mantissa) as an integer of `f`
* `e` power-of-two exponent of `f`
* `p` is the actual number of bits in `m`. This can be less than requested since any number of the low bits can be zero.
* `m` is adjusted by the number of trailing zeroes
* `e` is likewise adjusted

\\
and if you copied that into a file and ran it (`sollya whateverFileName`) you'd get:

	[|3.1415926535897966687116422690451145172119140625, 221069929750889, -46, 0|]
    [|3.1415926535897932384626433832794812270636909610957, 63719069007931157819013617823235, -104, 1|]
    [|0.61803398874989312616889947094023227691650390625, 173961102589770, -48, 1|]

\\
Taking the last one we have:

$$
\frac{2}{\sqrt{5}+1} \approx 0.61803... = 173961102589770 \cdot 2^{-47}
$$

\\
Finally we have the information we need to attack the factoring. If we take the integer and plug it into WolframAlpha it gives the prime factorization:

$$
173961102589770 = 2 \cdot 3 \cdot 5 \cdot 103 \cdot 25913 \cdot 2172581
$$

\\
We can set aside any powers of two for the moment (this is our trailing zeroes) and want to take the remaining terms and attempt to partition them into two sets where the product of each is small enough to fit in our floating point format without rounding. For floats that means they need to be $\le 2^{24}-1$. For this example we have two choices:

$$ \begin{align*}
 A & = 3 \cdot 2172581  \\ 
 B & = 5 \cdot 103 \cdot 25913 \\
 \\ 
 A & = 5 \cdot 2172581 \\
 B & = 3 \cdot 103 \cdot 25913 
\end{align*} $$

\\
Since we were successful with this integer we could have used the simpler sollya routine `add_k_q` to get the same results. There's no reason I can see to prefer one over the other. Toy example code:

{% highlight c %}
float add_inverse_phi(float x)
{
  const float A = ((3.f*2172581.f)/70368744177664.f);
  const float B = (5.f*103.f*25913.f);
  return fmaf(A,B,x);
}
{% endhighlight %}


\\
Now let's look at $\pi$ for *binary32* (an example usage is a part of range reduction for one and two parameter arc-tangents). The initial integer isn't going to work out this time so we're going to look at a range of value. Our target integer rounded up so we want to inspect: $ \left\\{I,I-1,I+1,I-2,I+2...\right\\} $. If it had rounded down we'd flip the sign on the ordering.

{: .center }
|     | $2^t$  | prime factorization         |
|:--- | :---:  | ---:                        |
| I   | 0 |$ 5867 \cdot 44279 \cdot 850973 $                 |
| I-1 | 3 |$ 317 \cdot \color{red}{87172685233}  $                        | 
| I+1 | 1 |$ 3 \cdot 5 \cdot 11 \cdot 1061 \cdot \color{red}{631393853} $ | 
| I-2 | 0 |$ 3 \cdot 7 \cdot 19 \cdot \color{red}{554059974313} $  | 
| I+2 | 0 |$ 13 \cdot 61 \cdot 73 \cdot 14879 \cdot 256661 $       |
| I-3 | 1 |$ 17 \cdot 4253 \cdot \color{red}{1528816543} $         |
| I+3 | 2 |$ \color{red}{55267482437723} $                         |


\\
The $2^t$ column is the extra powers of two and the red integers are too large, so we have two candidates in this window $I$ and $I+2$. I already spoiled that $I$ doesn't work out and that's because the smallest times the largest is too big. But we're good with $I+2$:

$$ \begin{align*}
 A & =  61 \cdot 256661 \\ 
 B & =  13 \cdot 73 \cdot 14879 \\
 \\
 \pi & \approx AB \cdot 2^{-46} \\
 \\
 \frac{\text{RN}\left( \pi \right) - \pi}{\pi} & \approx 2.78275 \times 10^{-8} & \approx 1.110 \times 2^{-26} \\
 \\
 \frac{AB \cdot 2^{-46} - \pi}{\pi} & \approx 1.01388 \times 10^{-14} & \approx 1.011 \times 2^{-47} 
\end{align*} $$

\\
And some toy code.

{% highlight c %}

static const f32_pair_t f32_mk_pi = {.h = (float)(61*256661), .l= (float)(13*73*14879)*0x1.0p-46f};

float whatever(float x)
{
  const float pi_a = f32_mk_pi.h;
  const float pi_b = f32_mk_pi.l;

  // do stuff like negate pi_b and/or maybe multiply it by 1/2,1/4...whatever
  // and probably some other stuff as well.

  return fmaf(pi_a, pi_b, r);
}
{% endhighlight %}

