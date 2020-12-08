# Dendritic Freezing Using Finite Elements
This project aims to model freezing by using Finite Element Meshes. Specifically we aim to model some common branching/dendritic phenomena we observe in winter sceneries.



![enter image description here](https://i.pinimg.com/originals/c8/ef/3b/c8ef3bab645b5d451d3e44e89b2c8322.jpg =300x300)

![enter image description here](https://i.imgur.com/XKAA6Wy.png =400x250)



 Freezing occurs when a water particle collides with the surface of a pre-existing ice structure, has the correct temperature profile to freeze (mainly that it is below the freezing temperature of water) *and* it can find an appropriate bonding spot to attach itself to said surface.  Unpacking this, we will see there are three main main contributions to the freezing process, the transportation and motion of a free water particle through space, the transportation of heat through space, and finally the ability for the particle to find a bonding location in the ice. 

Usually, freezing occurs in scenarios where one of these three processes are slower than the others. For example, if particles attach immediately to the surface with no difficulty and we are in a medium that conducts heat instantly, the entire process is limited by the particle motion in that medium. If we are in a super saturated medium (so particles are readily available everywhere) and heat conducts instantly through the medium, the freezing process is limited by the kinetics of the freezing process. 

As a result, we can either have *diffusion-limited* freezing or *kinetics-limited*. Diffusion-Limited freezing corresponds to freezing limited by one the first two processes (heat conduction and particle motion can both be modelled with diffusion). Diffusion processes are well studied and can lead to visually interesting branching growth patterns.  Kinetics-limited freezing corresponds to the freezing limited by the last process, where *kinetics* refers to the chemistry term for the rate at which a chemical or physical reaction occurs. Kinetics-Limited growth is responsible for some key features in some freezing processes, such as the flat and sharp faceted characteristics of snowflakes. Kinetics-Limited freezing is still widely debated in the physics literature today and there is no widely agreed upon mathematical model for it. As a result, our work will focus on the former Diffusion-Limited freezing, a process which is implicitly responsible for the unstable dendritic growth characteristic of frost.
## The Physics Equations
### Heat Diffusion
There are three principal ways in which heat transfer can occur: conduction, convection, radiation. We will assume convection and radiation are negligible. In real life, they usually aren't but we will show we can still reproduce interesting patterns with just diffusion. 
The traditional heat diffusion equation, as can be derived from conservation of energy, goes like this:
$$
 \frac{\partial T(\boldsymbol{x})}{\partial t} = C \Delta T(\boldsymbol{x}) + q(\boldsymbol{x}) \quad \quad \quad  \boldsymbol{x} \in \Omega
$$

Where $T(\boldsymbol{x})$ is the temperature at a point $\boldsymbol{x}$ in our domain $\Omega$, $C$ is a constant that includes specific heat capacity, density, and conductance constants, a negative sign, the $\Delta$ operator represents the Laplacian at that point, and $q(\boldsymbol{x})$ represents any sources or sinks present in our domain (maybe a point in space continuously emits heat, like a lightbulb). The equation above is a standard diffusion equation, and to gain intuition for it, consider a single (non-source) point in space. The temperature (or whatever quantity is diffusing) at that point in space will only change if the temperature at the points around it is different than its own. In this case,  according to Fourier's law,  $\boldsymbol{\phi} = - k \nabla T$, the heat flux $\phi$ points  from high temperatures to low temperatures, and the speed at which heat will go from high temperatures to low temperatures is modulated by the conductance $k$.  Maybe better intuition for this is if you think of the Laplacian as taking an average of neighboring locations, then the temperature at any point will change to match the average temperature at it's neighborhood. 
Furthermore, for the time being, let's just assume sources/sinks $q(\boldsymbol{x})$ are nonexistent, and instead we can model them by setting dirichlet boundary conditions to the above equation.

![enter image description here](https://i.imgur.com/3AMdJOq.gif =400x400)

### Interface Motion
We can think of our domain $\Omega$ as being split into two sub-domains $\Omega_s$ and $\Omega_l$ for the solid and liquid domains respectively. We denote $\Gamma$ to represent the interface of the two domains, where one domain becomes the other. For each of these sub-domains, we need to solve the diffusion equation in the previous section. However, we also need to figure out how our interface $\Gamma$ is moving.  This has been studied for literally hundreds of years, and the equation governing the motion of the interface is as follows:
$$
\boldsymbol{v}(\boldsymbol{x}) \cdot \boldsymbol{n}(\boldsymbol{x}) =D \frac{\partial T} {\partial \boldsymbol n(\boldsymbol{x})} =D (\frac{\partial T} {\partial \boldsymbol n(\boldsymbol{x})} |_s - \frac{\partial T} {\partial \boldsymbol n(\boldsymbol{x})} |_l)  \quad \quad \quad \boldsymbol{x} \in \Gamma
$$

Where $\boldsymbol{v}(\boldsymbol{x})$ represents the velocity for a point on our interface $\Gamma$, $\boldsymbol{n}$ is the normal at that same point on our interface (pointing from the solid domain to the liquid), D is a constant that includes conductance and latent heat of fusion. The above equation is often called the Stefan Condition and can once again be derived using conservation of energy. The stefan condition can also be understood intuitively: Consider the case where there is no phase change along the interface. The absence of phase change is another way of saying that $\boldsymbol{v}(\boldsymbol{x}) \cdot \boldsymbol{n} = 0$.  This implies that the heat flux going into the interface from the solid side $\frac{\partial T} {\partial \boldsymbol n}|_s$ is equal to the heat flux leaving the interface from the liquid side $\frac{\partial T} {\partial \boldsymbol n}|_l$. If the two are not equal to each other, then their difference must go towards either melting or freezing. That's all this equation is saying, the difference in heat fluxes along the interface is proportional to the speed at which the interface moves. The equation above is usually rewritten like below, to emphasize we only care about the magnitude of the speed along the normal direction, and that we can calculate the flux of temperature along the normal to be the temperature gradient, dotted with the normal.
$$
V_n=D (\nabla T |_s - \nabla T|_l)\cdot \boldsymbol{n(\boldsymbol{x})}  \quad \quad \quad \boldsymbol{x} \in \Gamma
$$
![enter image description here](https://i.imgur.com/SXMqDqG.png)
### Boundary Conditions and Melting/Freezing Temperature
To solve the coupled system of equations given by heat diffusion and the stefan condition, we need to impose boundary conditions on our diffusion equation for points along our interface. In other words, what should the temperature at the interface always be equal to. Obviously, it makes sense to set this to the bulk melting/freezing temperature of ice:
$$
T(\boldsymbol{x}) = T_M \quad \quad \quad \boldsymbol{x}  \in \Gamma
$$

Where $T_M = 273\degree K = 0\degree C$ for ice. This is a completely reasonable answer, and works well for modelling large scale freezing and melting processes, such as an ice sculpture melting, or the  freezing of a lake starting from the top and going downwards. However, many of the interesting visual frost phenomena we observe occur at a smaller scale, where surface tension of the ice plays a non-negligible role. You can think of surface tension as being a force present in the ice that tries to keep the surface as smooth as possible. Surface tension is easily modelled by the Gibbs-Thomson equation, which simply replaces the temperature equation above with the following:
$$
T(\boldsymbol{x}) = T_M(1 - \sigma \kappa )\quad \quad \quad \boldsymbol{x}  \in \Gamma 
$$

Where $\sigma$ represents the strength of the surface tension force and is often called the capillary length, $\kappa$ represents curvature (positive if the center of curvature is inside the solid, negative outside). Let's gain some intuition about this formula. Assume $\sigma$ is constant $\sigma = 0.002$. Imagine a 2D scenario where your interface is represented as a straight line, $\Gamma(\boldsymbol{x}) = 0$. Now let's say the interface is perturbed by a sine wave, so now your interface is given by the equation $\Gamma(\boldsymbol{x}) = sin(f\boldsymbol{x})$ where f represents frequency. You can convince yourself that $f$ implicitly encodes the curvature along the sine wave. 
 If you have a very low frequency $f \approx 0$, then you can see that the curvature along the peaks of your sin wave will be very very small. This means the temperature according to Gibbs-Thomson along your sine wave is close to the bulk melting temperature $T_M$. In fact, if the curvature is zero, we get back our planar interface where the curvature everywhere is zero and the melting temperature is then *exactly* the bulk melting temperature $T_M$. However, if you have a high frequency $f >> 0$, you'll have a lot of waves in a small area and your curvature will be very, very high along the peaks of your sin waves. Plugging this information into the Gibbs-Thomson equation, we can see that surface tension will have a non-negligible effect, ie the temperature at the interface will not be the same as the bulk melting temperature $T_M$, but will in fact be lower at the peaks, and higher at the crests.

### Quasi-Steady State Assumption
Before we talk about the discretization, we will simplify the diffusion equation in the first section:
$$
 \frac{\partial T(\boldsymbol{x})}{\partial t} = C \Delta T(\boldsymbol{x})  \quad \quad \quad  \boldsymbol{x} \in \Omega
$$

Because it takes much more energy to freeze a water molecule than to heat it up, or alternatively because in general the interface motion is MUCH slower than the heat diffusion process, we can assume that as soon as the ice melts, the heat released from that melt is immediately dissipated perfectly into it's environment. This allows us to reduce the diffusion equation above to the Laplace equation:
$$
\Delta T(\boldsymbol{x}) = 0 \quad \quad \quad  \boldsymbol{x} \in \Omega
$$

This is an assumption that is used very often in crystal growth physics/math literature, both for numerical solutions and analytical ones. This is a powerful assumption, as we can now say that the temperature field is harmonic in the domain. This now means, to determine the temperature at any point in time, we *only* need the boundary conditions of our problem, we no longer require the temperature field at previous points in time.

## The Discretization
In summary, these are the three equations that will govern our freezing process:
$$
\begin{cases}
1 \quad \Delta T= 0& \quad \quad \quad  \boldsymbol{x}  \in \Omega_s, \Omega_l \\
2 \quad V_n =D \frac{\partial T} {\partial \boldsymbol{n} }& \quad \quad \quad \boldsymbol{x} \in \Gamma \\
3 \quad T = T_M(1 - \sigma \kappa )&\quad \quad \quad \boldsymbol{x}  \in \Omega_s, \Omega_l
\end{cases}
$$

We will only be working in a 2D domain. We discretize our volumetric domain with triangles, where the interface is a set of edges shared by triangles in the solid and liquid domain. We assume our temperature field is piecewise linear over the 2D domain, a temperature is assigned to each triangle vertex(including vertices along the interface), and is interpolated using barycentric coordinate shape functions for non-vertex points in our domain.

![enter image description here](https://i.imgur.com/sbzPiHG.png)
Our overall algorithm for the whole simulation looks like this:

    For Each Timestep in the Simulation, do:
	    Solve Poisson Equation (Eq 1 and 3)
	    Move Interface In Time (Eq 2)
	    Remesh Domain 
	


### Solving Laplace Equation
We need to solve the Laplace equation on our triangle mesh, with appropriate boundary conditions given by the Gibbs Thomson equation, as well as any sources we wish to model in our environment.  We will assume that  our ice is fully submerged in a bounding box. We will set the temperature at each vertex on the boundary of our bounding box to have a set value of $-10 \degree C$. This is like imagining our little box is constantly pressed between 4 super cold solid plates whose temperature is fixed at -10. Then, for each vertex on the sold-liquid interface we will find the mean curvature. Note that the contour of our solid is the interface, and if our domain is discretized with a triangle mesh, then the contour is discretized with a polygonal curve. Calculating the curvature in 2D is trivial: to find the signed curvature at a given node, we can find the exterior angle made by it's two incident edges, as shown in assignment 7. 
![Calculating 2D discrete curvature, by using the exterior angle of a discrete curve](https://i.imgur.com/htcdp1m.jpg =400x200)

To solve the Laplace equation with these Dirichlet Boundary conditions, we will instead work with the equivalent problem of minimizing an energy functional. The corresponding energy functional that must be minimized is a little bit like an antiderivative to the Laplace equation:
$$
\boldsymbol v^* = \min_{\boldsymbol{v}} \boldsymbol v^T \boldsymbol C \boldsymbol v
$$

Where $\boldsymbol C$ is the cotan Laplacian matrix.  See the libigl tutorial on how to derive this! We can use `mqwf` from lilbigl to solve this energy minimization problem.

### Move Interface In Time
Moving the interface in time is a little more involved than the simple Laplace equation solve. Recall the motion of our mesh is determined by the following equation:
$$
V_n=D (\nabla T |_s - \nabla T|_l)\cdot \boldsymbol{n}  \quad \quad \quad \boldsymbol{x} \in \Gamma
$$

Where $D$ is a known constant. The main unknown we have in this equation is finding the gradient of our temperature field along every vertex in the interface. Note that because our shape functions used in our FEM discretization are piecewise linear, the gradient of the temperature field is *constant* along each triangle face. Furthermore, each edge in our interface has *two* triangle faces associated with it; one triangle face belonging to the solid domain, and one belonging to the liquid. Therefore the gradient quantities we are interested in the above equation, $\nabla T |_s$ and $\nabla T|_l$ are well defined along an interface edge. We calculate the gradient at each face using the $\boldsymbol G$ operator as defined by libigl and multiplying the temperature field by it. This will give a gradient vector for each face. So, we carefully have to go through each countouring edge in our interface, find the incident faces on that edge, index the temperature gradient vector, and dot the difference in gradients with the normal vector $\boldsymbol n$ associated with each edge. 

Note that this gives us a scalar $V_n$ quantity associated with each edge, indicating the velocity of each "point" along that edge would move along it's corresponding normal direction. This causes a problem for us however, because we can't directly move the edges, if we wish to maintain mesh connectivity,  we can only move the vertices! So we have a desired motion associated with each edge, and we wish to move our vertices to best match that motion. This is also known as a face offsetting problem in Geometry Processing. These can be solved naively in 2D by moving each edge independently by the offset, then recalculating it's intersection point with its neighboring edges. We take a different approach and instead formulate an energy we wish to minimize for the optimal vertex motion. The continuous version of our energy is given by :
$$
\min_{\overline{v}_i \forall i \in |V|} \int_\Gamma (\overline{v}_t(s)\cdot \hat{n} - \overline{v}(s) \cdot \hat{n})^2 ds
$$

Unpacking this, we want that for each point in our smooth *not-yet-discrete* contour to have a solution velocity $\overline{v}(s)$ that, when projected to the normal, be as close as possible to a target velocity $\overline{v}_t(s)$, where $s$ represents the parameter in our arclength parameterized curve. 
Discretizing  our curve, we says that  $\overline{v}(s)$ is piecewise linear along our curve, with velocities defined on vertices, while $\overline{v}_t(s)$ is piecewise constant along the curve, with each edge having a constant target velocity.
This results in the following discretized energy:
$$
\boldsymbol{V}^* =\min_{\boldsymbol{V}} \boldsymbol{V}^T\boldsymbol{M}\boldsymbol{V}-2\boldsymbol{V}_{pn}^T\boldsymbol{A}\boldsymbol{V} 
$$

Where $\boldsymbol{V}$ is a $(|N_v| \times 3)$ matrix, where each row represents a velocity vector we wish to solve for all vertices $N_v$. $V_{pn}$ is a $(|N_e| \times 1)$ vector that represents the target velocity projected to the normal at that point, in our case it would be the velocity given at each edge by our Stefan Condition, dotted with the normal of that edge, for all edges $N_e$. M is a $(|N_v| \times |N_v|)$ mass matrix derived in this document, and is calculated by integrating the shape functions along each edge, dotted with the normal. A is a $(|N_e| \times |N_v|)$ matrix that projects the influence of each edge onto the vertices to which it is incident, dotted with some normals. Detailed explanations on how to assemble these matrices is given in this pdf.
This energy is well suited for our problem because it is well defined along our discretized interface. Unlike other face offsetting problems, this energy keeps corners sharp and gives a well defined velocity for each vertex.

Finally, now that we have a velocity associated with each interface vertex, all that's left is for us to move it forward in time, with the forward update:
$$
\boldsymbol x^{i+1} = V_n^i \boldsymbol n^i + \boldsymbol x^i
$$
Where $i$ denotes the current timestep.

### Remeshing
Finally, after all that hard work, we note that if the interface is moved to much in time, some vertices will eventually move past others and the mesh will become "tangled". To avoid this, we take the simple approach and simply remesh by using the triangle remesher. The way this works is we retain the boundary polygon of our bounding box, as well as the polygon representing our interface. We retain the connectivity of these boundaries, and throw away all the other info regarding the interior of our domain. Then, inputting our old boundaries/interfaces into triangle as constraints, we tell triangle to remesh the entire interior given these boundaries. This ensures that we maintain a nice high quality mesh for every step in our simulation. Note also that we can throw away the interior of our domain specifically because we made the quasi-steady state assumption. If we kept the standard heat equation, we would need to retain the interior as, to find the temperature at one node, we'd need to know it's temperature from the previous timestep as well. Because of our quasi-steady state assumption, we say that the energy from the previous timestep is fully dissipated at the next timestep, and therefore to recalculate it we only need the boundary temperatures.

This also implies that we don't even need to *mesh* our domain in the first place. If the temperature at each point only depends on the boundary, why spend all that time remeshing. If you were thinking this, you'd be absolutely right, and this is the next step in this project, to make this work with the boundary element method. Note that a similar problem to this (with the same equations) has already been solved with BEM in this paper. 

![Suggested boundary only discretization](https://i.imgur.com/yOFe26L.png)
### Final Notes
You might be thinking that nothing we have modelled here explains why dendritic growth happens. The answer is that you're right, nothing mentioned here explicitly models unstable dendritic growth, but rather this unstable growht arises implicitly from these equations. It turns out that when solving the heat equation, the temperature gradient is higher in areas with high curvature! Because interface motion is directly proportional to the temperature gradient, regions of high curvature will therefore move faster than areas of low curvature. This only makes our geometry have even HIGHER curvature the next timestep, and the process repeats. That is how this instability occurs. 

| ![enter image description here](https://media.giphy.com/media/5McygFmeVdQk1Bd3Mi/giphy.gif =600x300) 
![enter image description here](https://media.giphy.com/media/KgZqEFOfiL3qCjzozQ/giphy.gif)
## How to Run The Assignment Code
All the code is written in matlab. You need gptoolbox as well as the triangle extension. To find details on how to install the triangle extension, check it out here. If you're on windows, note I attached a precompiled triangle binary to my submission. Just go to `path_to_triangle.m` and change the path defined in variable `s` to be where `s = ./`
 
 To run the code, simply run either `dendriticGrowth_box.m` , for an ice crystal submerged in a liquid box as described above, or `dendriticGrowth_line.m` for a slightly different and modified problem to the one above, but most of the physics is still the same.


Note you can vary around various parameters at the top of each of these files to change the conditions. These parameters are surface tension, timestep length, bounding box temperature (set this to 10 if you want to see melting), maximum triangle area, and some more. Note that some configurations will crash our remesher.

