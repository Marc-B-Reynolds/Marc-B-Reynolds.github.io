---
layout:       post
title:        Quaternion are really Complex
tagline:      
categories:   [quaternions]
tags:         [foundation]
description:  some quaternion wankery support material for other posts
---

There are a number of quaternion topics I'm thinking of writing about which might be easier to understand, depending on your math background, if we take some non-standard viewpoints.  The purpose of these post is to sketch out some of these for later reference.  Forgive my use and massive abuse of notation.  The title is part pun, part self-mockery.  The only reason to read anything here is for a slightly more detailed hand-waving version of my hand-waving somewhere else.  The less detailed hand-waving version elsewhere is intended to be the understandable version BTW.

<div class="alert alert-success" role="alert" markdown="1">
If you cannot follow the formalism or if it is making your eyes bleed too much, I will give plane (da dum dum tish!) English summaries marked in green.
</div>

$ \newcommand{\hcomplex}[1]{\mathbb{C}_{#1}} $
$ \newcommand{\pwrap}[1]{\left( #1 \right)} $
$ \newcommand{\set}[1]{\left\{ #1 \right\}} $

------

Preliminaries <small>the really boring part</small>
------

\\
The short version of this section is: $\forall Q \in \mathbb{H}$ the {1, Q}-plane forms a commutative subalgebra of $\mathbb{H}$ isomorphic to $\mathbb{C}$.


\\
Quaternions ($\mathbb{H}$) were formulated to be an extension of complex number ($\mathbb{C}$) from two spatial dimensions to three.  We can write a complex number as:

$$ \begin{eqnarray}
  z & = & a + b~\mathbf{i} \label{zsb} \\ 
    & = & \left(a,b\right) \label{zcoord}
\end{eqnarray} $$

where $a,b\in \mathbb{R}$ and $\mathbf{i}$ is a unit bivector.  We will call "$ a $" the scalar part and "$ b~\mathbf{i} $" the bivector part and associate $\mathbb{C}$ with a two dimensional real valued vector space ($\mathbb{R}^2$).  As such we will consider $z$ to be the coordinate: $ \eqref{zcoord} $.

Complex numbers require two axiomatic definitions:

* The product of the unit bivector: $\mathbf{i}^2=-1$
* The (complex) conjugate:  $ z^* = \left(a + b~\mathbf{i}\right)^* = a - b~\mathbf{i} $

\\
all other properties follow from these.


\\
Quaternions follow the same form as complex numbers and we can write a quaternion as a scalar/bivector pair:

$$ \begin{eqnarray}
  Q & = & a + b~\mathbf{u} \label{qsb} \\ 
    & = & \left[a,b\right] \label{qcoord}
\end{eqnarray} $$

where $a,b\in \mathbb{R}$ and $\mathbf{u}$ is a unit bivector, however a bivector in three dimensions requires three scalars instead of one to be represented. This brings the total number to four so we associate $\mathbb{H}$ with four dimensional real valued vector space ($\mathbb{R}^4$).  However given a known/implied unit bivector we can also denote a quaternion as the order pair of reals: $ \eqref{qcoord} $.

Quaternions require a single axiomatic definition: the product of two bivectors. The quaternion conjugate can be formed with products and additions.  Without defining either they work out to be:

* The square of any unit bivector: $\mathbf{u}^2=-1$
* The conjugate:  $ Q^* = \left(a + b~\mathbf{u}\right)^* = a - b~\mathbf{u} $

\\
So informally the subset of quaternions with fixed $\mathbf{u}$ is identical (isomorphic) to $\mathbb{C}$.  We can associate the set of all unit bivectors of $\mathbb{H}$ with the three dimensional unit sphere:

$$ \begin{equation} \label{S2}
\mathbb{S}^{2}=\left\{ \mathbf{u}\in\mathbb{R}^{3}:\left|\mathbf{u}\right|=1\right\} 
\end{equation} $$

\\
and the subset of quaternions with fixed $\mathbf{u}$ can be denoted as:

$$ \begin{equation} \label{Cu}
\hcomplex{u}=\left\{ \left(a + b~\mathbf{u}\right) \in\mathbb{H}:a,b~\mathbb{\in R}~; \mathbf{u}\in\mathbb{S}^{2}\right\} 
\end{equation} $$

\\
The sets $ \hcomplex{u} $ and $ {\mathbb{C}\_{-u}} $ are logically distinct since they have opposite orientations, but are they are equivalent in that they represent the same set of coordinates (oriented vs. unoriented planes).  The full set of quaternions is (equivalent to) the set of all these complex planes:

$$
\mathbb{H}=\left\{ X \in \hcomplex{u} : \mathbf{u} \in\mathbb{S}^{2} \right\} 
$$

\\
Informally allowing the unit bivector to vary every quaternion can be expressed as $ \eqref{qsb} $.

<div class="alert alert-success" role="alert" markdown="1">
There is a plane that spans the line of reals and any non-real ($b\neq 0$) quaternion. Collectively the quaternions in this plane behave like the standard set of complex numbers.  Every non-real quaternion falls in exactly one (or two depending on your chosen perspective) of these complex planes.  All of these complex plane intersect at the line of reals.
</div>

\\
A quaternion can also be denoted as a coordinate in ${\mathbb{R}^4}$ :

$$ \begin{equation} \label{qxform}
Q = a + \left(x,~y,~z\right) 
\end{equation} $$

\\
As practitioners we are working with quaternions over reals (modeled over floating-point) in (typically) form $ \eqref{qxform} $ and not quaternions over arbitrary expressions that evaluate to reals.  Therefore our form $ \eqref{qsb} $ always have $ b = \left| \pwrap{x,~y,~z} \right| \geq 0 $.  We are always working with the positive half-plane:

$$ \begin{equation} \label{CuPlus}
\hcomplex{u}^{+} =\left\{ \left(a + b~\mathbf{u}\right) \in\hcomplex{u} : b\geq 0\right\} 
\end{equation} $$

<br>

* $Q$ is always in the positive half plane of the {1,Q}-plane (the implied $\hcomplex{u}$).
* All rotations about $u$ are in $\hcomplex{u}$ (origin excluded).
* The projected line $sQ$, where $s$ is a non-zero scalar, all represent the same rotation.
* Unit quaternions in $\hcomplex{u}$ form a unit circle in the plane.
* The projected line $sQ$ intersects the unit-sphere (unit-circle in the plane) at $Q$ and $-Q$, if $Q$ is a unit quaternion. Negating a quaternion negates the implied angle and the orientation of the unit bivector.

<br>

### NOTES:

Although there are some conceptual presentation and notational differences in the above there is only one contradiction with norm.  Specifically I am defining a quaternion to be scalar/bivector pair instead of a scalar/vector pair.  Apparently Simon L. Altmann provides a proof in *"Rotations, Quaternions and Double Groups"*. An informal argument is as follows:

* $\mathbb{C}$ and $\mathbb{H}$ cannot be sub/supersets.  A 2D vector requires two reals so $\mathbb{C}$ is not a scalar/vector pair (that would be three reals).
* Crazy physicist call 2D bivectors pseudoscalars (1 real) and 3D bivectors pseudovectors or axial-vectors (3 reals).
* Both $\mathbb{C}$ and $\mathbb{H}$ naturally describe angular and not directional quantities.  Angular quantities such as angular velocity, acceleration, etc. are bivectors (unless you're a crazy physicist).
* Converting a 2/3D bivector to linear alegbra yields skew-symmetric (antisymmetric) matrices, which is consistent the homo/isomorphisms of both to linear algebra.
* Bivector notion is consistent with Clifford/Geometric algebra.

Also note that the complex conjugate is an involution: $ {(z_0 z_1)^{\*}={z_0}^{\*}{z_1}^{\*}} $ and the quaternion conjugate $ {(Q_0 Q_1)^{\*}={Q_1}^{\*}{Q_0}^{\*}} $ is an anti-involution.  However if $Q_0$ and $Q_1$ are in the same complex plane, then product commutes and they behave the same in the plane.

<br>

------

Quaternion valued complex functions <small>bear with me</small>
------
{:#complexFunction}

\\
If we take the set of all definable functions over $\mathbb{C}$ which only use real valued coefficients, then over quaternions they are elements of:

$$ {\mathbf{\mathbb{C}}_{u}}^{\mathbf{\mathbb{C}}_{u}} $$

\\
which we will term "*(quaternion valued) complex maps or functions*".  So a function from this set $f$ is the map:

$$ f:\mathbf{\mathbb{C}}_{u}\rightarrow\mathbf{\mathbb{C}}_{u} $$

\\
In simpler terms if $f$ with complex input yields:

$$ f\left(a + b~\mathbf{i}\right) = c + d~\mathbf{i}$$

then with quaternion input is:

$$ f\left(a + b~\mathbf{u}\right) = c + d~\mathbf{u}$$

Some trival examples are addition/subtraction/multplication by reals and the conjugate.  More interesting is all real analytic functions (all smooth functions? don't know..doesn't matter).  Notably exp, log, real powers and (well) 2D conformal maps.

There is an applied gotcha here when the input is real valued ($b=0$) and the result is complex valued ($d\neq 0$).  In the complex case there is only one plane and the result is single valued. In the quaternion case we either have a known $\mathbf{u}$ in which it is single valued in that plane, otherwise we have an infinite number of results in every complex plane.  As an example, the complex result:

$$\sqrt{-1} = \sqrt{\left(-1,0\right)} = \left(0,1\right) $$

in quaternions requires a known plane, otherwise the result is the set of all unit bivectors (the unit sphere $\mathbb{S}^{2}$).

As an aside the Mandelbrot and Julia set formulations are *"complex maps"* and this is why directly translating to quaternion versions are doomed to failure.  The result is the same in every complex plane of $\mathbb{H}$ so the results are a surface-of-revolution about the line of reals.

Alone this is of minimal (at best) interest, but combined with other properties of quaternions can become interesting. As a non-random example if we want reproduce the behavior that some function $f$ has in the standard complex-plane to the plane spanned by {0,A,B}, (where the real value is mapped to $A$):

$$ f\pwrap{BA^{*}}~A = A~f\pwrap{A^{*}B} $$


<div class="alert alert-success" role="alert" markdown="1">
There are a bunch of functions over quaternions that map an input in some complex plane to the same complex plane.  These functions are effectively 2D.
</div>



<br>

------

Moving out of the plane
------

\\
The product of two quaternions \$ {A \in \hcomplex{u}} \$ and  \${B \in \hcomplex{v}} \$ can be expanded as:

$$ \begin{eqnarray*}
  AB & = & \left(a + b~\mathbf{u}\right) \left(c + d~\mathbf{v}\right)   \\ 
     & = & ac + bc~\mathbf{u} + ad~\mathbf{v} + bd~\mathbf{uv}
\end{eqnarray*} $$

\\
To complete we need the rule for the product of two unit vectors $\mathbf{uv}$.  The axiomatic rules for the product are to be found everywhere and the bivector product works out to be:

$$ 
\mathbf{uv} = \mathbf{u} \times \mathbf{v} - \mathbf{u}\cdot\mathbf{v}
$$

\\
Completing the expansion and collecting the scalar and bivector parts yields:

$$ 
AB = \left(ac - bd \pwrap{\mathbf{u}\cdot\mathbf{v}} \right) + \left( bc~\mathbf{u} + ad~\mathbf{v} + bd \pwrap{\mathbf{u} \times \mathbf{v}}  \right) 
$$

<div class="alert alert-danger" role="alert" markdown="1">
This space left intentionlly blank.  I ran out of steam for the moment.
</div>


