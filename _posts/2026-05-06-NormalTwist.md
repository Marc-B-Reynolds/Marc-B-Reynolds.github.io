---
layout:       post
title:        'A normal/twist transform for TBN representation'
categories:   [math]
tags:         [graphics, tangent]
description:  'Covers a transform that converts the tangent and bitangent of a 3D frame to a point on the unit circle.'
plotly:       true
---

\\
Yeah...I know. Exactly what the world needs: *"yet another representation of a rotation"*.  So I'll spitball it first and show the math later.

Given a normal, tangent and bitangent (which I'll associate with Z,X,Y respectively) transform the tangent & bitangent information to a point on a unit circle.  I'll assume the use-case is quantization which needs to encode/decode a 3D normal and 2D normal (lots of off-the-shelf options for both). This factorization allows one to choose the accuracy of $\mathbf{n}$ and $\mathbf{t}$ & $\mathbf{b}$ (as a pair) independently. Obviously if we sometimes only need $\mathbf{n}$ then there's no additional work since it's untransformed. The forward and inverse transforms each require a divide (sad-face), the forward uses two scalar products and the inverse uses 8 to produce both $\mathbf{t}$ & $\mathbf{b}$.  Recall that transforming a quaternion representation back to TBN requires 12 scalar products and no divides (using the most commonly seen method). 

Other than showing the code and breaking down the math...that's it.

NOTE: I'm not claiming this is novel.  It's inconceivable that it is but I've not run across a write-up. The quaternion based derivation might be though.

<br>

------

Limited input version
------

\\
Let's start by looking at a simplified version of the two transforms. These versions are good when $n_z \ge 0$. As $n_z$ becomes increasingly negative the errors will grow until it completely explodes.  Stable versions are in the next section.


```c?encoding=ascii
// unstable version: will blow up as n.z is approaching -1.
// for n.z on [0,1] produces identical results to the stable version.
// ever increasing error as z moves more negative.

vec2f_t normal_twist_limited_encode(vec3f_t n, vec3f_t t)
{
  float s = t.z/(n.z + 1.f);
  float x = fmaf(n.x,s,-t.x);
  float y = fmaf(n.y,s,-t.y);

  return vec2f(x,y);
}

static inline void normal_twist_limited_decode(vec3f_t* t, vec3f_t* b, vec3f_t n, vec2f_t e)
{
  // spitball: 1 divide and 8 scalar products
  float s  =  n.y/(n.z+1.f);
  float tz =  fmaf(e.x,n.x, e.y*n.y);   // xx+yy
  float bz =  fmaf(e.x,n.y,-e.y*n.x);   // xy-yx
  float sx =  s*tz;                     // s(xx+yy)
  float sy =  s*bz;                     // s(xy-yx)
  float tx = -fmaf(e.x,n.z,sy);         // s(yx-xy)-xz
  float bx =  fmaf(e.y,n.z,sx);         // s(xx+yy)+yz
  float ty =  sx-e.y;                   // s(xx+yy)-e.y
  float by =  sy-e.x;                   // s(xy-yx)-e.x

  t[0] = vec3f(tx,ty,tz);
  b[0] = vec3f(bx,by,bz);
}
```

<br>

------

Stable version, defs and helper funcs
------

\\
Extending to stable computation for all normals requires sprinkling around some sign flips.

The code style here is to show the flow of computations and not be concise nor optimized for
any specific architecture. This listing can be found at [godbolt](https://godbolt.org/z/9jErq1dGa). 

```c?encoding=ascii

// clang specific for example code to look (kinda) shader like
// and to have values passed in register. The GCC extension
// uses array access notion which makes it harder to
// follow as example code IMHO.
typedef float vec2f_t __attribute__((ext_vector_type(2)));
typedef float vec3f_t __attribute__((ext_vector_type(4)));


static inline vec2f_t vec2f(float  x, float  y)            { return (vec2f_t){x,y}; }
static inline vec3f_t vec3f(float  x, float  y, float  z)  { return (vec3f_t){x,y,z,0}; }


static inline uint32_t f32_to_bits(float x)
{
  uint32_t u; memcpy(&u, &x, 4); return u;
}

static inline float f32_from_bits(uint32_t x)
{
  float f; memcpy(&f, &x, 4); return f;
}

// isolate sign bit of 'a'
static inline uint32_t f32_sign_bit(float a)
{
  return f32_to_bits(a) & 0x80000000;
}

// if 'v' is float and 's' is all clear (except sign bit)
// return copysign(1.f,s)*v
static inline float f32_mulsign(float v, uint32_t s)
{
  return f32_from_bits(f32_to_bits(v)^s);
}


// Given normal 'n' transform tangent 't' to a unit circle.
// The transform is logically: find the swing rotation that
// rotates +/-Z (sign is that of z component of 'n') to 'n' and
// use that to rotate 't' to the XY plane.  The result is
// equivalent to an implied twist.
//
// The selecting of +/-Z makes the transform numerically
// stable.

vec2f_t normal_twist_encode(vec3f_t n, vec3f_t t)
{
  uint32_t sb = f32_sign_bit(n[2]);     // isolate sign bit of n.z
  float    sz = f32_mulsign(1.f,sb);    // sgn(n.z) = copysignf(1.f,n.z)

  float s  = t.z/(n.z + sz);
  float x  = fmaf(n.x,s,-t.x);
  float y  = fmaf(n.y,s,-t.y);

  x = f32_mulsign(x,sb);
  y = f32_mulsign(y,sb);

  return vec2f(x,y);
}


// Given normal 'n' and encoded twist 'e': produce 't' and 'b'
void normal_twist_decode(vec3f_t* t, vec3f_t* b, vec3f_t n, vec2f_t e)
{
  // spitball: 1 divide and 8 scalar products to apply
  // the twist to 'n' to restore t & b

  uint32_t sb = f32_sign_bit(n[2]);     // isolate sign bit of n.z
  float    sz = f32_mulsign(1.f,sb);    // sgn(n.z) = copysignf(1.f,n.z)
  float    ey = f32_mulsign(e.y,sb);
  float    ex = f32_mulsign(e.x,sb);

  float s  =  n.y/(n.z + sz);
  float tz =  fmaf(e.x,n.x, e.y*n.y);
  float bz =  fmaf(e.x,n.y,-e.y*n.x);
  float sx =  s*tz;
  float sy =  s*bz;
  float tx = -fmaf(e.x,n.z,sy);
  float bx =  fmaf(e.y,n.z,sx);
  float ty =  sx-ey;
  float by =  sy-ex;

  // bitangent needs to be negated if (n.z < 0)
  bx = f32_mulsign(bx,sb);
  by = f32_mulsign(by,sb);
  bz = f32_mulsign(bz,sb);
  
  t[0] = vec3f(tx,ty,tz);
  b[0] = vec3f(bx,by,bz);
}


// Given normal 'n' and encoded twist 'e': produce 't'
//   mainly to show that only producing one of t or b is about
//   the same cost as producing both.
vec3f_t normal_twist_decode_t(vec3f_t n, vec2f_t e)
{
  // spitball: only removes 1 scalar product to directly generate one
  uint32_t sb = f32_sign_bit(n[2]);
  float    sz = f32_mulsign(1.f,sb);
  float    ey = f32_mulsign(e.y,sb);

  float s  =  n.y/(n.z + sz);
  float tz =  fmaf(e.x,n.x, e.y*n.y);
  float bz =  fmaf(e.x,n.y,-e.y*n.x);
  float sx =  s*tz;
  float sy =  s*bz;
  float tx = -fmaf(e.x,n.z,sy);
  float ty =  sx-ey;
  
  return vec3f(tx,ty,tz);
}
```

<br>

------

Math
------

\\
The math is a simple extension of previous blog posts. Start with ["Orthonormal basis from normal via quaternion similarity"](https://marc-b-reynolds.github.io/quaternions/2016/07/06/Orthonormal.html) which generates a $\mathbf{t}$ and $\mathbf{b}$ from $\mathbf{n}$ we have expression(s) for a *torque minimal* rotation ($\pm~\mathbf{z}$ to $\mathbf{n} $):

<div style="font-size: 120%;">
$$
Q_s = \left\{ \begin{array}{llll}
    \frac{1+n_z}{\sqrt{2+2n_z}} + \left(\frac{-n_y}{\sqrt{2+2n_z}},~\frac{n_x}{\sqrt{2+2n_z}},~0\right) && n_z \ge 0   \\
    \frac{1-n_z}{\sqrt{2-2n_z}} + \left(\frac{n_y}{\sqrt{2-2n_z}}, ~\frac{-n_x}{\sqrt{2-2n_z}},~0\right) && n_z \lt 0  \\
    \frac{1+\text{sgn}(n_z)~n_z}{\sqrt{2+2~\text{sgn}(n_z)~n_z}} + \left(\frac{-\text{sgn}(n_z)~n_y}{\sqrt{2+2~\text{sgn}(n_z)~n_z}},~\frac{\text{sgn}(n_z)~n_x}{\sqrt{2+2~\text{sgn}(n_z)~n_z}},~0\right) && \text{generalized}
\end{array}\right.
$$
</div>

\\
We can consider this rotation (about an axis in the $\mathbf{xy}$ plane) to be the *swing* portion of a [swing twist decomposition](https://marc-b-reynolds.github.io/quaternions/2022/01/31/QuatAxisFactor.html). Rotating the tangent by the conjugate of $Q_s$ (opposite angle) removes the *swing* producing a point on the unit circle:

$$
\mathbf{t}' = Q^{*}_{s}~\mathbf{t}~Q_s = \left(t_x-\frac{n_x t_z}{n_z+\text{sgn}(n_z)},~t_y-\frac{n_y t_z}{n_z+\text{sgn}(n_z)},~0\right)
$$

\\
Translating `normal_twist_encode` back into an expression gives:

$$
\mathbf{e} = \text{sgn}(n_z)\left(\frac{n_x t_z}{n_z+\text{sgn}(n_z)}-t_x,~\frac{n_y t_z}{n_z+\text{sgn}(n_z)}-t_y\right)
$$

\\
Notice that it's the expression $\mathbf{t}'$ above which has been negated and then negated again if $n_z \lt 0$. These modifications
are using advanced knowledge of how the reverse expressions will turn out.  To perform the reverse transform we simply need to
rotate back in place:


$$
\begin{array}{llll}
\mathbf{t} & = Q_s~\left(\phantom{-}e_x,e_y,0\right)~Q^{*}_{s} \\
\mathbf{b} & = Q_s~\left(-e_y,e_x,0\right)~Q^{*}_{s}
\end{array}
$$

\\
I think it easier to look at the expansions of the two cases and merge the differences so expanding for $n_z \ge 0$

$$
\begin{array}{llll}
\mathbf{t} & = \left(\frac{n_y}{n_z+1}\left(e_y n_x-e_x n_y\right)-e_x n_z,
                    ~\frac{n_y}{n_z+1}\left(e_x n_x+e_y n_y\right)-e_y,
                    ~e_x n_x+e_y n_y\right) & & n_z \ge 0 \\
\mathbf{b} & = \left(\frac{n_y}{n_z+1} \left(e_x n_x+e_y n_y\right)+e_y n_z,
                    ~\frac{n_y}{n_z+1} \left(e_x n_y-e_y n_x\right)-e_x,
                    ~e_x n_y-e_y n_x\right) & & n_z \ge 0
\end{array}
$$

\\
And likewise for $n_z \lt 0$ case. I'm writting an expression for $-\mathbf{b}$ because our unrotated right handed frame is $ \set{-\mathbf{z},\mathbf{x},-\mathbf{y}} $ but showing the expression for the left handed $ \set{-\mathbf{z},\mathbf{x},\mathbf{y}} $ is easier to eye-ball the differences.

$$
\begin{array}{rlll}
\mathbf{t} &  = \left(
                     \frac{n_y}{n_z-1}\left(e_y n_x-e_x n_y\right)-e_x n_z,
                    ~\frac{n_y}{n_z-1}\left(e_x n_x+e_y n_y\right)+e_y,
                    ~e_x n_x+e_y n_y\right) & & n_z \lt 0 \\
-\mathbf{b} & = \left(
                     \frac{n_y}{n_z-1} \left(e_x n_x+e_y n_y\right)+e_y n_z,
                    ~\frac{n_y}{n_z-1} \left(e_x n_y-e_y n_x\right)+e_x,
                    ~e_x n_y-e_y n_x\right) & & n_z \lt 0
\end{array}
$$

\\
And finally merge the two cases:

$$
\begin{array}{llll}
\mathbf{t} &  = \phantom{\text{sgn}(n_z)}\left(
                     \frac{n_y}{n_z+\text{sgn}(n_z)}\left(e_y n_x-e_x n_y\right)-e_x n_z,
                    ~\frac{n_y}{n_z+\text{sgn}(n_z)}\left(e_x n_x+e_y n_y\right)-\text{sgn}(n_z)~e_y,
                    ~e_x n_x+e_y n_y\right) \\
\mathbf{b}  & = \text{sgn}(n_z)\left(
                     \frac{n_y}{n_z+\text{sgn}(n_z)} \left(e_x n_x+e_y n_y\right)+e_y n_z,
                    ~\frac{n_y}{n_z+\text{sgn}(n_z)} \left(e_x n_y-e_y n_x\right)-\text{sgn}(n_z)~e_x,
                    ~e_x n_y-e_y n_x\right)
\end{array}
$$
