---
layout:       post
title:        Table of basic reversible integer operations
categories:   [math]
tags:         [integer]
description:  A table and some comments on integer bijections.
---

This is a small table of reversible integer functions that can be implemented in a small number of opcodes. This is adapted from a similar table by Bret Mulvey[^pscarab] which I ran across from a Pixar technical note[^pixar].  I'm basically going to fill in what I think are some interesting points.  I'm going to assume all operations are on whole registers, otherwise you need to sprinkle in AND masks as needed.

## Basic table

The table is intended to rather minimal.  Some forward notes:

* Mostly avoiding duplicate entries so no need to list the inverse functions seperately.
* The *field* indicates the mathematical field where the operation is most easily expressed.  I'm using $\mathbb{Z}$ to indicate integers modulo some power of two[^integers] and $\mathbb{F_2}$ for informally bit-wise operations.

<br>

{: .center }
|     field     | $f\left(x\right)$ |$f^{-1}\left(x\right)$|notes| 
| :---:         |  ---              | --- |
|$\mathbb{Z} $  | `x  = -x`         | `x  = -x`                    |   |
|$\mathbb{F}_2$ | `x ^= k`          | `x ^= k`                     |   |
|$\mathbb{Z} $  | `x += k`          | `x -= k`                     |   |
|$\mathbb{Z} $  | `x *= k`          | `x *= mod_inverse(k)`        | 1, k must be odd |
|$\mathbb{Z} $  | `x += x << k`     | `x *= mod_inverse(1+(1<<k))` | 1,2 |
|$\mathbb{Z} $  | `x -= x << k`     | `x *= mod_inverse(1-(1<<k))` | 1,2 |
|$\mathbb{F}_2$ | `x ^= x >> k`     |  --- see below               | 3 |
|$\mathbb{F}_2$ | `x ^= x << k`     |  --- see below               | 3,5 |
|$\mathbb{F}_2$ | `x = rotl(x,k)`   | `x = rotr(x,k)`              | 4 |
|$\mathbb{F}_2$ | `x = permute(x)`  |  --- see below               | 3,4 |
|$\mathbb{F}_2$ | ---  various      |  --- see below               | 6   |

<br>

1. Where `mod_inverse` is the multiplicative inverse[^modinverse].
2. Ok the forward operations are multiplication.
3. See below for inverse function.
4. Bit rotations are a special case of bit permutations. Other examples include byte-swaps and bit reversals.
5. Special case of "carryless multiplication".
6. See below.

Since all of these operations are invertible[^bijection], any composition of them is likewise invertible.  So if we perform three in a row:

$$ \begin{eqnarray*}
x' &=& f_2\left(f_1\left(f_0\left(x\right)\right)\right) \\
   &=& \left(f_2 \circ f_1 \circ f_0\right)\left(x\right)
\end{eqnarray*} $$

\\
then the inverse is simply *solving* the equation which is chaining the inverse functions in the opposite order.


$$ \begin{eqnarray*}
x &=& f^{-1}_0\left(f^{-1}_1\left(f^{-1}_2\left(x'\right)\right)\right) \\
   &=& \left(f^{-1}_0 \circ f^{-1}_1 \circ f^{-1}_2\right)\left(x'\right) \\
   &=& \left(f_2 \circ f_1 \circ f_0\right)^{-1}\left(x'\right)
\end{eqnarray*} $$

\\
Numbers aside: The total number of invertible functions is equal to the number of permutations.  For a bit-width of $\beta$ there are $2^{\beta}!$ unique invertible functions.  A big number but most are boring and worse not fast to compute.

<br>

## Quick bit-wise operations as $\mathbb{F_2}$

\\
We can consider a bit to be an integer modulo 2 and can describe basic bit operations as the two element finite field[^simp] $\mathbb{F}_2$. As a field we need to support the operations: addition, subtraction, multiplication and division.  If have two computer integer values `a` and `b` where the lowest bit of each are the values we're performing operations on:


{: .center }
|   $\mathbb{F}_2$  operation  | C-like operation | C-like reduced |notes| 
| ---     |  ---                | ---   | 
| $-a$      |  `(-a) & 1`       | `a`   | since $1+1=0$, $-1=1$ also $-0=0$   |
| $a+b$     |  `(a+b) & 1 `     | `a^b` |   |
| $a-b$     |  `(a-b) & 1 `     | `a^b` | or via negate: $a+\left(-b\right) = a+b$  |
| $a~b$     |  `(a*b) & 1 `     | `a&b` |   |
| $a^{-1}$  |  `mod_inverse(a)` | `a`   |   |

\\
In turn a collection of bits (like a register) can treated as a vector of $\mathbb{F}_2$ elements.  So the *C-like reduced* are the vector operations (multiplication and multiplicative inverse taken to be component-wise).  We can go a step further and consider our vectors as column matrices which can be transformed by a matrix.  Taking 32-bit registers then we can represent a transform by a 32x32 matrix.  Well we can't represent vector addition this way, we'd have to increase the matrix size to 33x33, ignoring this here since XORing by a constant isn't super interesting.  Taking the convention that our 32-bit integer `x` in matrix form is $x=\left(x_0,x_1,\ldots\right)^{T}$ where $x_0$ is the low bit, then there are two basic transform matrices we need to consider: permutations and shifts (bit rotations can be represented by either).  Transforming $x$ to $x'$ becomes:


$$ x' = M~x$$

\\
And since all the functions are reversible then the inverse of $M$ exists.  And from group theory (Cayley's theorem[^cayley]) we have:

$$M^n = I$$

\\
where $n$ is some positive integer.  The smallest $n$ is the (finite) order (period/cardinality) of $M$. For computation we don't necessarily need the smallest $n$ if the op in question is $I$ for any greater power.  This gives us:

$$ M^{-1}~M = M^n \implies M^{-1} = M^{n-1} $$


\\
Aside: If the dimensions of our matrix is $\beta$ (using example of $\beta=32$) then there are $2^{\beta^2}$ matrices (so example is $2^{1024}$).  The number of invertible matrices is: 

$$\prod _{i=0}^{n-1} \left(2^\beta-2^i\right) = 2^{\beta^2}\text{QPochhammer}\left(2^{-\beta}, ~2, ~\beta\right)$$

\\
The percentage of invertible matrices approaches a limit of 28.8788%. (Mathematica told me about the RHS above.  QPochhammer[^qsymbol] is news to me.)


<br>

### Bit permutations

Bit permutations can be represented by, well, a permutation matrix[^permmatrix] which has exactly one '1' entry in each column and row and everything else is zero. Permutation matrices are orthogonal so if $P$ is a permutation matrix then $P^{-1} = P^{T}$.  The identity matrix $I$ qualifies.  Bit-twiddling resources are loaded with examples, like the public available *"Matters Computational"*[^arndt].  Since I mentioned numbers before, we have $\beta!$ of these (big, but very small compared to total number of invertible).  Some examples:

* Identity matrix
<p align="center"><canvas id="ident" width="161" height="161"></canvas></p>
* Exchange matrix[^exchange] $J$ represents bitreversal. It's a self inverse so $J^2=I$.
<p align="center"><canvas id="brev" width="161" height="161"></canvas></p>
* Various other *swap* operations, notably *byteswap* (shown below), which are also involutions.
<p align="center"><canvas id="bswap" width="161" height="161"></canvas></p>
* Rotating left can be represented by a circulant matrix[^circulant], specifically the cyclic-permutation matrix let's call that $C$. Raising $C$ to a positive integer $k$ is a left rotate by $k$, negative $k$ becomes right rotate and $C^{\beta} = I$ (any $k \bmod \beta = 0$ as well, showing order). 
<p align="center"><canvas id="rot" width="161" height="161"></canvas></p>
<p align="center"><code class="highlighter-rouge" id="rott">PDEP</code></p>

* BMI2 ops parallel bits deposit (`PDEP`) and parallel bits extract (`PEXT`) can be used to build bit permutations.


<br>

### Bit shifts

Given the stated conventions then a right shift of one bit is represented by an upper shift matrix[^shiftm] let's call that $R$ and a left shift of one bit by lower shift matrix $L$. These alone are not invertible since we have information loss of the shifted away bits.  Raise each to positive integer power `k` is a shift by `k` bit positions and if $k \geq \beta$ the result becomes zero (your computer language may not like this in code). These two matrices are mutual transposes.


### Xorshifts

The right xorshift: `x = x ^ (x>>k)` which RHS translates into $x + R^k~x = \left(I + R^k\right)x $, so the matrix is:

$$ M = \left(I + R^k\right) $$

\\
if $k\geq \beta$ then $M=I$.

<p align="center"><canvas id="xsr" width="161" height="161"></canvas></p>
<p align="center"><code class="highlighter-rouge" id="xsrt">PDEP</code></p>

\\
The left xorshift: `x = x ^ (x<<k)` likewise becomes the transformation matrix:

$$ M = \left(I + L^{k}\right) $$

<p align="center"><canvas id="xsl" width="161" height="161"></canvas></p>
<p align="center"><code class="highlighter-rouge" id="xslt">PDEP</code></p>


### Inverse of xorshifts
{:#xorshiftInverse}

\\
If you don't happen to know how to find the inverse of an xorshift then it's pretty straightword to understand. Let's just look at some code for the standard *gray code*:

{% highlight c %}
uint32_t u32_to_gray(uint32_t x) { return x^(x >> 1); }

uint32_t gray_to_u32(uint32_t x)
{
  x = x ^ (x >>  1);
  x = x ^ (x >>  2);
  x = x ^ (x >>  4);
  x = x ^ (x >>  8);
  x = x ^ (x >> 16);
  return x;
}
{% endhighlight %}

\\
After the *binary to gray* conversion the top bit is left unchanged. To reverse we repeat the forward transform and we're back to two bits, next step we've got four, etc.  Verifying all of the bits get restored is tedious but not very difficult.  For any right xorshift by `k`, the inverse is a sequence of $i$ right xorshifts with shifts of: $2^i~k$ (each twice of the previous...we're double number of restored bits with each), where $i$ is enough that the next shift would be greater-than or equal to the bit width.  Examples for three, five and nine:

{% highlight c %}
// extra ugly names to (hopefully) fit on screen
uint32_t irxorshift_3(uint32_t x) { x=x^(x>>3); x=x^(x>> 6); x=x^(x>>12); return x^(x>>24); }  // next would be 48
uint32_t irxorshift_5(uint32_t x) { x=x^(x>>5); x=x^(x>>10); return x^(x>>20); }               // next would be 40
uint32_t irxorshift_9(uint32_t x) { x=x^(x>>9); return x^(x>>18); }                            // next would be 36
{% endhighlight %}

\\
If the shift amount is half the width or greater then the xorshift is an involution.

Let's play the matrices a bit. Mechanically expand a composition of two right xorshifts `x ^= x>>a; x^= x>>b` :


$$ \begin{eqnarray*}
 \left(I + R^b\right)\left(I + R^a\right) 
   & = & I^2 + I R^a + I R^b + R^aR^b \\
   & = & I + R^a + R^b + R^{a+b}       \\
   & = & \left(I + R^a\right)\left(I + R^b\right) 
\end{eqnarray*} $$

\\
so same direction xorshifts commute. The wikipedia list of gray to binary runs in the opposite order of that above. For dead horse thing this expression is the same as above: `x = x ^ (x>>a) ^ (x>>b) ^ (x>>(a+b))`

Setting $a=b=k$ and plugging in gives squaring:

$$ \begin{eqnarray*}
 \left(I + R^k\right)^2 & = & I + 2 R^k + R^{2k} \\
                        & = & I + R^{2k}

\end{eqnarray*} $$

which give the 'bit-doubling' operation we used above for the inverse.  Combine this with the Cayley's theorem result and we can directly form the inverse.


 $$ M^{-1} = M^{n-1} = \left(I + R^k\right)\left(I + R^{2k}\right)\left(I + R^{4k}\right)\ldots $$

\\
I go into the details a bit more in the companion XorRotate post[^xorrot].  Everything said about right xorshift apply to left.

<br>

### Carryless multiplies and their inverses

\\
Carryless multiplies opcodes are now[^lowlatency] fast enough to be usable on Intel-a-likes.  Let's just jump in with a quick foolishly naive software implementation for 32x32 mod 2^32:


{% highlight c %}
// Don't do this at home: intended for clarity
// PCLMUL* opcodes
uint32_t cl_mul_32(uint32_t a, uint32_t b)
{
  uint32_t r = 0;

  while (b != 0) {
    if (b & 1) { r ^= a; } a <<= 1; b >>= 1;
  }

  return r;
}
{% endhighlight %}

\\
If we were to change `r ^= a` to `r += a` we'd have a really bad integer product flipping to long multiplication in modulo $\mathbb{Z}$ from modulo $\mathbb{F}_2$ (where addition is XOR).  If we were to fix either to a constant (let's choose `b`) then we could expand this based on the bits of `b` which are set: `(x<<b0)^(x<<b1)...`  where `bn` is the position of the n<sup>th</sup> set bit.  You might be shocked to discover that if the lowest bit is set then this carryless product has a multiplicative inverse.  The forward transform in matrix would then look like:


$$
M = I + L^{b_1} + L^{b_2} + \ldots 
$$

\\
our now old-friend the left xorshift. And to find a generic inverse we simply need to compute $M^{\beta-1}$.


{% highlight c %}

// slightly less painful forward version
uint32_t cl_mul_32(uint32_t a, uint32_t b)
{
  __m128i v1 = _mm_cvtsi64_si128(a);
  __m128i v2 = _mm_cvtsi64_si128(b);
  __m128i r  = _mm_clmulepi64_si128(v1,v2, 0);
  return (uint32_t)_mm_cvtsi128_si64(r);
}

// modulo multiplicative inverse of carryless multiplier 'x'
// cl_mul_32(k, cl_mul_inv_32(k)) = 1 for all odd 'x'.
uint32_t cl_mul_inv_32(uint32_t x)
{
  uint32_t r = x;             // M

  // generate and compose: M^2,M^4,M^8,M^{16}
  for (int i=0; i<4; i++) {
    x = cl_mul_32(x,x);      // square previous power of M
    r = cl_mul_32(r,x);      // compose
  }

  // returning M^{-1} = M^{31}
  return r; 
}
{% endhighlight %}

<br>

The code above is structured for clarity and will produce unneeded operations. An actually implementation for current Intel intrinics would look like:

{% highlight c %}

uint64_t cl_mul_inv_32(uint32_t v)
{
  __m128i x = _mm_cvtsi32_si128(v);
  __m128i r = x;

  for(int i=0; i<4; i++) {
    x = _mm_clmulepi64_si128(x,x, 0);
    r = _mm_clmulepi64_si128(r,x, 0);
  }
  return _mm_cvtsi128_si64(r);
}

uint64_t cl_mul_inv_64(uint64_t v)
{
  __m128i x = _mm_cvtsi64_si128(v);
  __m128i r = x;

  for(int i=0; i<5; i++) {
    x = _mm_clmulepi64_si128(x,x, 0);
    r = _mm_clmulepi64_si128(r,x, 0);
  }
  return _mm_cvtsi128_si64(r);
}

{% endhighlight %}


<br>

### Others

* Basics of XOR rotates are in another post[^xorrot].
* AND XOR rotates (not motived yet)
* CRC32-C (not motived yet)

<br>

------

References and Footnotes
------

[^integers]:   Mentally replace that by $\mathbb{Z}_{2^b}$ or $\mathbb{Z}/2^{b}$ if you're offended.
[^modinverse]: *"Integer multiplicative inverse via Newton's method*", 2017 ([page]({{site.base}}/math/2017/09/18/ModInverse.html))
[^pscarab]:    *"Hash functions"*, Bret Mulvey ([page](http://papa.bretmulvey.com/post/124027987928/hash-functions))
[^pixar]:      *"Correlated Multi-Jittered Sampling"*, Andrew Kensler 2013 ([PDF](http://graphics.pixar.com/library/MultiJitteredSampling/paper.pdf))
[^permmatrix]: Wikipedia: Permutation matrix ([page](http://en.wikipedia.org/wiki/Permutation_matrix))
[^simp]:       Well $\mathbb{F}_{2}$ is more general than math on bits...don't need any of that here.
[^shiftm]:     Wikipedia: Shift Matrix ([page](http://en.wikipedia.org/wiki/Shift_matrix))
[^bijection]:  Wikipedia: Bijection ([page](http://en.wikipedia.org/wiki/Bijection))
[^exchange]:   Wikipedia: Exchange matrix ([page](http://en.wikipedia.org/wiki/Exchange_matrix))
[^cayley]:     Wikipedia: Cayley's theorem ([page](http://en.wikipedia.org/wiki/Cayley%27s_theorem))
[^qsymbol]:    Wikipedia: Q-Pochhammer ([page](http://en.wikipedia.org/wiki/Q-Pochhammer_symbol))
[^arndt]:      *"Matters Computational"*, Jorg Arndt, 2010 ([page](http://www.jjj.de/fxt/))
[^circulant]:  Wikipedia: Circulant matrix ([page](https://en.wikipedia.org/wiki/Circulant_matrix#Properties))
[^xorrot]:     *Basic XOR-rotates and their inverse*, 2017 ([page]({{site.base}}/math/2017/10/13/XorRotate.html))
[^lowlatency]: Carryless multiplies on Intel Haswell+ and AMD Jaguar+ microarchitecture. SEE: Fog's [tables](http://www.agner.org/optimize/instruction_tables.pdf)

<script>

// can't be bothered
function popcount(x) {var m=1,r=0; while (m>0) {if(x&m){r++;} m = (m<<1)&0xffffffff;} return r;}

function grid_draw_d(ctx,x,y,len)
{
  x = 5*x+1;
  y = 5*y+1;
  for (var i=0; i<len; i++) {ctx.fillRect(x, y, 4, 4);x+=5;y+=5;}
}

function grid_draw_j(ctx)
{
  var x = 31*5+1, y=1;
  for (var i=0; i<32; i++) { ctx.fillRect(x, y, 4, 4); x-=5;y+=5;}
}

function grid_draw_bswap(ctx)
{
  grid_draw_d(ctx, 24,  0, 8);
  grid_draw_d(ctx, 16,  8, 8);
  grid_draw_d(ctx,  8, 16, 8);
  grid_draw_d(ctx,  0, 24, 8);
}

function grid_draw_ld(ctx,pos) { grid_draw_d(ctx,0,pos,32-pos); }
function grid_draw_ud(ctx,pos) { grid_draw_d(ctx,pos,0,32-pos); }

function grid_draw_rot(grid)
{
  var ctx = grid.ctx;
  var pos = grid.k;
  grid_draw_ld(ctx,pos);
  grid_draw_ud(ctx,32-pos);

  document.getElementById("rott").textContent = "rot(x," + pos + ")";
  grid.k = (grid.k+1)&0x1f;
}

function grid_draw_xsr(grid)
{
  var ctx = grid.ctx;
  var pos = grid.k;
  grid_draw_ud(ctx, 0);
  grid_draw_ud(ctx, pos);

  document.getElementById("xsrt").textContent = "x^(x>>" + pos + ")";
  grid.k++; if (grid.k==32) {grid.k=1;}
}

function grid_draw_xsl(grid)
{
  var ctx = grid.ctx;
  var pos = grid.k;
  grid_draw_ud(ctx, 0);
  grid_draw_ld(ctx, pos);

  document.getElementById("xslt").textContent = "x^(x<<" + pos + ")";
  grid.k++; if (grid.k==32) {grid.k=1;}
}

function grid_clear(ctx,size)
{
  var e = (size-1)/5;
  //ctx.fillStyle = 'rgb(230,240,255)';
  //ctx.fillRect(0, 0, size, size);
  ctx.fillStyle = 'rgb(255,255,255)';
  ctx.fillRect(0, 0, size, size);
  ctx.fillStyle = 'rgb(230,240,255)';
  for (var j=0; j<e; j++) {
    for (var i=0; i<e; i++)  {
      ctx.fillRect(5*i+1, 5*j+1, 4, 4);
	}
  }
}

function grid_draw(grid)
{
  var ctx  = grid.ctx;
  var size = 5*grid.dim+1;

  grid_clear(ctx,161);
  //ctx.fillStyle = grid.bgfill;
  //ctx.fillRect(0, 0, size, size);
  
  ctx.fillStyle = grid.fgfill;
  grid.fdraw(grid);

  window.setTimeout(grid.draw, 200);
}

function grid_make(name,func)
{
  var c = document.getElementById(name);
  var x = c.getContext('2d');

  var grid = {
	ctx:    x,
	draw:   function() { grid_draw(grid); },
	fdraw:  func,
	k:      6,
	dim:    32,
	cnt:    0,
	bgfill: 'rgb(230,240,255)',
	fgfill: '#f00'
  };

  return grid;
}

function grid_draw_static(name,func)
{
  var id  = document.getElementById(name);
  var ctx = id.getContext('2d');
  grid_clear(ctx,161);
  ctx.fillStyle = '#f00';
  if (func != null) func(ctx);
}

grid_make('rot', grid_draw_rot).draw();
grid_make('xsr', grid_draw_xsr).draw();
grid_make('xsl', grid_draw_xsl).draw();

grid_draw_static('ident', function(ctx) { grid_draw_ld(ctx,0); });
grid_draw_static('brev',  function(ctx) { grid_draw_j(ctx); });
grid_draw_static('bswap', function(ctx) { grid_draw_bswap(ctx); });

</script>







