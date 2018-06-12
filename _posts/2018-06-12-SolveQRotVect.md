---
layout:       post
title:        Solving for quaternion rotation from a vector viewpoint
tagline:      
categories:   [quaternions]
tags:         [rotations, graphics]
description:  A short description of solving a quaternion rotation using vectors
---

Finding the quaternion that rotates unit vectors $\mathbf{a}$ to $\mathbf{b}$ is an operation I'm using over-and-over on this blog and it's really silly that I never explictly described it from standard vector viewpoint.  Just spewing out the standard geometeric defs of the cross and dot products:

$$
\begin{eqnarray*}
\mathbf{a} \times \mathbf{b} &=&  \lVert \mathbf{a} \rVert~ \lVert \mathbf{b} \rVert~ \sin\left(\theta\right)~ \mathbf{n} \\
\mathbf{a} \cdot \mathbf{b} &=&  \lVert \mathbf{a} \rVert~ \lVert \mathbf{b} \rVert~ \cos\left(\theta\right)
\end{eqnarray*}
$$

\\
where $\lVert \ . \rVert~ $ denotes magnitude and since we're limited ourselves to unit vectors reduces to:

$$
\begin{eqnarray*}
\mathbf{a} \times \mathbf{b} &=& \sin\left(\theta\right)~ \mathbf{n} \\
\mathbf{a} \cdot \mathbf{b} &=&  \cos\left(\theta\right)
\end{eqnarray*}
$$

\\
with $\mathbf{n}$ the unit vector orthogonal to the plane that spans $\mathbf{a}$ and $\mathbf{b}$ with right-handed orientation from $\mathbf{a}$ to $\mathbf{b}$.  Now if think of $\mathbf{n}$ as the axis of rotation and $\theta$ as the angle in an axis-angle representation then the rotation is the minimum magnitude angle that performs this operation (all the others can be described as a post composition of a rotation about $\mathbf{b}$). 

To create a quaternion we can simply look-up the often cited axis-angle to quaternion equation:

$$
Q = \cos\left(\frac{\theta}{2}\right) + \sin\left(\frac{\theta}{2}\right) \mathbf{n}
$$

\\
If it wasn't for that pesky "divide the angle by two" bit we'd be done.  Luckily we can just apply trig half-angle identities to make everything fit:

$$
Q = \frac{1+\mathbf{a} \cdot \mathbf{b}}{\sqrt{2+2~\mathbf{a} \cdot \mathbf{b}}} + \frac{\mathbf{a} \times \mathbf{b}}{\sqrt{2+2~\mathbf{a} \cdot \mathbf{b}}}
$$

\\
The quaternion equivalent of computing the cross and dot products and shoving the result into a quaternion can be performed by:

$$
\mathbf{b}\mathbf{a}^* = \cos\left(\theta\right) + \sin\left(\theta\right)~ \mathbf{n}
$$

\\
Since conjugation negates the angle and the product sums the resulting angle is: $\func{\text{angle}}{\mathbf{b}} - \func{\text{angle}}{\mathbf{a}}$

