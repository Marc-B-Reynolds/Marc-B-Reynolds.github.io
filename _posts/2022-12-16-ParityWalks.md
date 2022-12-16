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
// * increasing order odd  parity subset OEIS: A001969 (aka odious numbers) [2^31, 2^32)
// sign bit selects the parity, lower 31 is 'n' for the given a(n) 
// (P.S. see later)
uint32_t eop_code_32(uint32_t n)
{
  uint32_t r = bit_prefix_sum_32(n);   // cr_mul_32(r,-1) : inverse gray code
  return r ^ (r<<1);                   // cl_mul_32(r, 3) : bit-reflected gray code
}

// inverse function
uint32_t eop_code_inv_32(uint32_t n)
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

<br>

### Apparently it's easy to be evil <small>added 5 minutes later</small>

\\
What I should have done prior to pushing this post was to examine how the functions `eop_code_32` and `eop_code_inv_32` look in matrix form (same conventions as other [posts](https://marc-b-reynolds.github.io/math/2017/10/13/IntegerBijections.html)) respectively:

<p align="center" vertical-align="middle">
<canvas id="mat1"></canvas>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<canvas id="mat2"></canvas>
</p>

\\
Taking the matrix of `eop_code_32` (the left):

* the top row is all ones and "row times column" so the low bit of the output is the sum of all the input bits which is the parity
* all we have left is the first lower subdiagonal is all ones which corresponds to a left shift by one.

and the inverse function is pretty much the same song: top bit is the parity of the input and a left shift by one. Updated code:

{% highlight c %}

uint32_t eop_code_32(uint32_t n)
{
  return bit_parity_32(n) ^ (n<<1);
}

uint32_t eop_code_inv_32(uint32_t n)
{
  return (bit_parity_32(n)<<31) ^ (n>>1);
}

{% endhighlight %}

\\
This makes for a quite nice reduction: [godbolt](https://gcc.godbolt.org/z/1Ka7sjPda)

<script>

// colors for binary matrix visualization
//const mat_color = [[245,250,255],[50,0,0]];
const mat_color = [[230,240,255],[245,0,05]];

const cell_w = 4;  // drawn cell size (pixels) w/o padding

const mat1 = [
0xffffffff, 0x00000001, 0x00000002, 0x00000004, 0x00000008, 0x00000010, 0x00000020, 0x00000040, 0x00000080, 0x00000100, 0x00000200, 0x00000400, 0x00000800, 0x00001000, 0x00002000, 0x00004000, 0x00008000, 0x00010000, 0x00020000, 0x00040000, 0x00080000, 0x00100000, 0x00200000, 0x00400000, 0x00800000, 0x01000000, 0x02000000, 0x04000000, 0x08000000, 0x10000000, 0x20000000, 0x40000000];

const mat2 = [
0x00000002, 0x00000004, 0x00000008, 0x00000010, 0x00000020, 0x00000040, 0x00000080, 0x00000100, 0x00000200, 0x00000400, 0x00000800, 0x00001000, 0x00002000, 0x00004000, 0x00008000, 0x00010000, 0x00020000, 0x00040000, 0x00080000, 0x00100000, 0x00200000, 0x00400000, 0x00800000, 0x01000000, 0x02000000, 0x04000000, 0x08000000, 0x10000000, 0x20000000, 0x40000000, 0x80000000, 0xffffffff];

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


function viz_matrix_32(name,data)
{
  var w = 32;
  var h = 32;
  var cw  = (cell_w+1)*w+1;
  var ch  = (cell_w+1)*h+1;
  var id  = document.getElementById(name);
  var ctx = id.getContext('2d');
  var img = ctx.createImageData(cw,ch);
  var pix = img.data;
  var si  = 0;
  
  id.width  = cw;
  id.height = ch;

  pix.fill(255);

  for (var y=0; y<32; y++) {
    var bits = data[y];
    for (var x=0; x<32; x++)  {
	  var c = mat_color[bits & 1];
	  bits >>= 1;
	  viz_fill(pix, x, y, w, c[0],c[1],c[2]);
    }
  }

  ctx.putImageData(img,0,0)
}

viz_matrix_32('mat1', mat1);
viz_matrix_32('mat2', mat2);



