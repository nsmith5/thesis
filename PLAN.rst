

Goals for my PhD
================

This document is here to outline, in detail, the goals that I would like to accomplish during my PhD. Goals are split into major goals and the minor goals that make up the major goals. 


1.0 Demonstrate the Flexibility of the Binary XPFC model
=========================================================

The first major goal, and probably the simpliest is to demonstrate that with the new :math:`\epsilon(T)(c-c_0)^2/2` term in the free energy we should be able to produce all invariant reactions that a binary system can ever exhibit. The point here is not to dig into the intimate details of each of these reactions but to show, at a high level, that each one of them *can* be modelled.


1.1 The Syntectic Reaction :math:`l_1 + l_2 \rightarrow \alpha` 
----------------------------------------------------------------

To demonstrate the syntectic reaction we must complete the following tasks:

- Construct a free energy with the appropriate properties
- Construct a phase diagram
- Demonstrate that the system numerically. 
- Show nucleation of the :math:`\alpha` phase in simulation
- Research relevant systems and experiments to cite

Note: The spinodal in the liquid is *not* accurate. The coarsening of such a system requires hydrodynamics to be accurately described

1.2 The Monotectic Reaction :math:`l_1 \rightarrow l_1 + \alpha`
---------------------------------------------------------------- 

To demonstrate the monotectic reaction we must complete the following tasks:

- Construct a free energy with the appropriate properties
- Construct a phase diagram
- Demonstrate that the system numerically. 
- Research relevant systems and experiments to cite

1.3 The Eutectoid Reaction :math:`\gamma \rightarrow \alpha + \epsilon`
-----------------------------------------------------------------------

To demonstrate the eutectoid reaction we must complete the following tasks:

- Construct a free energy with the appropriate properties
- Construct a phase diagram
- Demonstrate that the system numerically.
- Research relevant systems and experiments

Note that there are different definitions of eutectoid in the literature

1.4 The Peritectoid Reaction :math:`\alpha + \epsilon \rightarrow \delta`
--------------------------------------------------------------------------

To demonstrate the eutectoid reaction we must complete the following tasks:

- Construct a free energy with the appropriate properties
- Construct a phase diagram
- Demonstrate that the system numerically.
- Research relevant systems and experiments

1.5 Solid State Spinodal :math:`\alpha \rightarrow \alpha_1 + \alpha_2`
------------------------------------------------------------------------

To demonstrate the eutectoid reaction we must complete the following tasks:

- Construct a free energy with the appropriate properties
- Construct a phase diagram
- Demonstrate that the system numerically.
- Research relevant systems and experiments

2.0 Explore the Eutectic Nucleation Pathway
===========================================

Here we want to look at specific kinetic pathway in a lot of detail. The pathway is the nucleation of a lamellar eutectic when a metastable solid solution :math:`\gamma^\prime` exists below the eutectic point. Simulations have shown that this state seems to be stabilized by elastic effects and we want to elucidate that mechanism in detail. We also want to show what processes lead to the dynamic decomposition of the :math:`\gamma^\prime` state in a growing front. 

2.1 Construct a Free Energy 
---------------------------

A concrete free energy needs to be developed with a phase diagram that is appropriate for the system. That is, we need a eutectic system in which a metastable solid solution is stabilized below the eutectic. The free energy must be careful in having a nucleation barrier in the amplitude not a spinodal. What is the effect of the T0 curve here 

2.2 Demostrate the Nucleation of :math:`\gamma^\prime`
------------------------------------------------------ 

Show that :math:`\gamma^\prime` all by itself at moderate undercoolings. Perhaps a TTT curve should be constructed to show the time scale of this type of nucleation. Don't get sidetracked trying to make theory for homogeneous nucleation. 

2.3 Demonstrate that a Growing Front can break down into Lamellar Eutectic
---------------------------------------------------------------------------

It is suspected that the mechanism for this is fluctuations in the concentrations at the interface. Is this true?

2.4 Show with an Amplitude model that Elastic Energy can Stabilize :math:`\gamma^\prime`
-----------------------------------------------------------------------------------------

Either by using the Elder model that has already been published or by deriving a new model from your theory, justify the fact that elastic energy is playing a role in the stability of the metastable state that can not be seen from the bulk free energy surface.

2.5 Show, if possible, that all of this can explain the Anomolous Eutectics of AgCu
------------------------------------------------------------------------------------

This amounts to looking at what happens under recalescense of the metastable and lamellar composite system. The metastable state is the least stable so it should melt first, leading to regrowth of the lamellar system back into the melt


3.0 Multiscale Phase field Crystal Model
=========================================

This work is completely seperate from the first two major goals. The foundations of a multiscale phase field crystal model have already been laid down, the basic idea is that the reference density should be promoted to a long wavelength field and both fields should be dynamical quantities. A supporting idea for this is that the correlation functions can now be considered as expansions of any long wavelength reference density!

3.1 Construct a Hard Sphere Inspired Model
-----------------------------------

Here are the tasks for the Hard Sphere Model

- Use the Percus-Yevick closure to make a direct correlation function
- Do a local density approximation of the hard sphere free energy for the reference free energy

3.2 Construct a Van der Waals Inspired Model
---------------------------------------------










