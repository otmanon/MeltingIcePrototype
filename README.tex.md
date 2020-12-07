# Freezing Using Finite Elements
This project aims to model freezing by using Finite Element Meshes.  Freezing occurs when a water particle collides with the surface of a pre-existing ice structure, has the correct temperature profile to freeze (mainly that it is below the freezing temperature of water) *and* it can find an appropriate bonding spot to attach itself to said surface.  Unpacking this, we will see there are three main main contributions to the freezing process, the transportation and motion of a free water particle through space, the transportation of heat through space, and finally the ability for the particle to find a bonding location in the ice. 

Usually, freezing occurs in scenarios where one of these three processes are slower than the others. For example, if particles attach immediately to the surface with no difficulty and we are in a medium that conducts heat instantly, the entire process is limited by the particle motion in that medium. If we are in a super saturated medium (so particles are readily available everywhere) and heat conducts instantly through the medium, the freezing process is limited by the kinetics of the freezing process. 

As a result, we can either have *diffusion-limited* freezing or *kinetics-limited*. Diffusion-Limited freezing corresponds to freezing limited by one the first two processes (heat conduction and particle motion can both be modelled with diffusion). Diffusion processes are well studied and can lead to visually interesting branching growth patterns.  Kinetics-limited freezing corresponds to the freezing limited by the last process, where *kinetics* refers to the chemistry term for the rate at which a chemical or physical reaction occurs. Kinetics-Limited growth is responsible for some key features in some freezing processes, such as the flat and sharp faceted characteristics of snowflakes. Kinetics-Limited freezing is still widely debated in the physics literature today and there is no widely agreed upon mathematical model for it. As a result, our work will focus on the former Diffusion-Limited freezing.
## The Physics Equations
### Heat Diffusion
There are three principal ways in which heat transfer can occur: conduction, convection, radiation. We will assume convection and radiation are negligible. In real life, they usually aren't but we will show we can still reproduce interesting patterns with just diffusion. 
The traditional heat diffusion equation, as can be derived from conservation of energy, goes like this:
$$
 \frac{\partial T(\boldsymbol{x})}{\partial t} = C \Delta T(\boldsymbol{x}) + q(\boldsymbol{x}) \quad \quad \quad  \boldsymbol{x} \in \Omega
$$

Where $T(\boldsymbol{x})$ is the temperature at a point $\boldsymbol{x}$ in our domain $\Omega$, $C$ is a constant that includes specific heat capacity, density, and conductance constants, a negative sign, the $\Delta$ operator represents the Laplacian at that point, and $q(\boldsymbol{x})$ represents any sources or sinks present in our domain (maybe a point in space continuously emits heat, like a lightbulb). The equation above is a standard diffusion equation, and to gain intuition for it, consider a single (non-source) point in space. The temperature (or whatever quantity is diffusing) at that point in space will only change if the temperature at the points around it is different than its own. In this case,  according to Fourier's law,  $\boldsymbol{\phi} = - k \nabla T$, the heat flux $\phi$ points  from high temperatures to low temperatures, and the speed at which heat will go from high temperatures to low temperatures is modulated by the conductance $k$.  Maybe better intuition for this is if you think of the Laplacian as taking an average of neighboring locations, then the temperature at any point will change to match the average temperature at it's neighborhood. 
Furthermore, for the time being, let's just assume sources/sinks $q(\boldsymbol{x})$ are nonexistent, and instead we can model them by setting dirichlet boundary conditions to the above equation.
### Interface Motion
We can think of our domain $\Omega$ as being split into two sub-domains $\Omega_s$ and $\Omega_l$ for the solid and liquid domains respectively. We denote $\Gamma$ to represent the interface of the two domains, where one domain becomes the other. For each of these sub-domains, we need to solve the diffusion equation in the previous section. However, we also need to figure out how our interface $\Gamma$ is moving.  This has been studied for literally hundreds of years, and the equation governing the motion of the interface is as follows:
$$
\boldsymbol{v}(\boldsymbol{x}) \cdot \boldsymbol{n}(\boldsymbol{x}) =D \frac{\partial T} {\partial \boldsymbol n(\boldsymbol{x})} =D (\frac{\partial T} {\partial \boldsymbol n(\boldsymbol{x})} |_s - \frac{\partial T} {\partial \boldsymbol n(\boldsymbol{x})} |_l)  \quad \quad \quad \boldsymbol{x} \in \Gamma
$$
Where $\boldsymbol{v}(\boldsymbol{x})$ represents the velocity for a point on our interface $\Gamma$, $\boldsymbol{n}$ is the normal at that same point on our interface (pointing from the solid domain to the liquid), D is a constant that includes conductance and latent heat of fusion. The above equation is often called the Stefan Condition and can once again be derived using conservation of energy. The stefan condition can also be understood intuitively: Consider the case where there is no phase change along the interface. The absence of phase change is another way of saying that $\boldsymbol{v}(\boldsymbol{x}) \cdot \boldsymbol{n} = 0$.  This implies that the heat flux going into the interface from the solid side $\frac{\partial T} {\partial \boldsymbol n}|_s$ is equal to the heat flux leaving the interface from the liquid side $\frac{\partial T} {\partial \boldsymbol n}|_l$. If the two are not equal to each other, then their difference must go towards either melting or freezing. That's all this equation is saying, the difference in heat fluxes along the interface is proportional to the speed at which the interface moves. The equation above is usually rewritten like below, to emphasize we only care about the magnitude of the speed along the normal direction, and that we can calculate the flux of temperature along the normal to be the temperature gradient, dotted with the normal.
$$
V_n=D (\nabla T |_s - \nabla T|_l)\cdot \boldsymbol{n(\boldsymbol{x})}  \quad \quad \quad \boldsymbol{x} \in \Gamma
$$

### Boundary Conditions and Melting/Freezing Temperature
To solve the coupled system of equations given by heat diffusion and the stefan condition, we need to impose boundary conditions on our diffusion equation for points along our interface. In other words, what should the temperature at the interface always be equal to. Obviously, it makes sense to set this to the bulk melting/freezing temperature of ice:
$$
T(\boldsymbol{x}) = T_M \quad \quad \quad \boldsymbol{x}  \in \Gamma
$$

Where $T_M = 273\degree K = 0\degree C$ for ice. This is a completely reasonable answer, and works well for modelling large scale freezing and melting processes, such as an ice sculpture melting, or the  freezing of a lake starting from the top and going downwards. However, many of the interesting visual frost phenomena we observe occur at a smaller scale, where surface tension of the ice plays a non-negligible role. You can think of surface tension as being a force present in the ice that tries to keep the surface as smooth as possible. Surface tension is easily modelled by the Gibbs-Thomson equation, which simply replaces the temperature equation above with the following:

$$
T(\boldsymbol{x}) = T_M(1 - \sigma \kappa )\quad \quad \quad \boldsymbol{x}  \in \Gamma $$

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
\Delta T= 0& \quad \quad \quad  \boldsymbol{x}  \in \Omega_s, \Omega_l \\
V_n =D \frac{\partial T} {\partial \boldsymbol{n} }& \quad \quad \quad \boldsymbol{x} \in \Gamma \\
T = T_M(1 - \sigma \kappa )&\quad \quad \quad \boldsymbol{x}  \in \Omega_s, \Omega_l
\end{cases}
$$

We will only be working in a 2D domain. We discretize our volumetric domain with triangles, where the interface is a set of edges shared by triangles in the solid and liquid domain. We assume our temperature field is piecewise linear over the 2D domain, a temperature is assigned to each triangle vertex(including vertices along the interface), and is interpolated using barycentric coordinate shape functions for general points in our domain.