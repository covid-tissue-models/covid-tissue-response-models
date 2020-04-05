
.. _title_start:

===============================
Model of Viral Tissue Infection
===============================

.. note::

    To cite this model please use the following:

    Josua Aponte-Serrano, T.J. Sego, Juliano F. Gianlupi, James A. Glazier,
    "Model of Viral Tissue Infection",
    https://github.com/covid-tissue-models/covid-tissue-response-models/tree/master/CC3D/Models/BiocIU/SARSCoV2MultiscaleVTM

.. _title_end:

.. _model_description_start:

Model Description
=================

This model describes select interactions between generalized epithelial
and immune cells and their extracellular environment associated with
viral infection and immune response at the cellular and intracellular
levels, and in the context of spatiotemporal dynamics. While at the
moment, the model is a generic viral infection model our hope is to
collaboratively develop it into a model of SARS-CoV-2 tissue infection
and Covid-19 progression. As such, it is intended to serve as a base
model for constructing and implementing more advanced models of targeted
cellular- and intracellular-level phenomena in tissue after initial
exposure. In its current state, it has not been formally peer-reviewed,
and should not be used for patient diagnostics or predicting clinical
outcomes. Rather, the model and its implementation can be used to
develop and interrogate mechanistic hypotheses about the spread of a
virus and how the interplay between viral spreading and immune response
determine the outcome of the disease, such as:

-  Why does the progression of the disease seem to be dependent on the
   initial viral exposure level?

-  Why is the start time of symptoms and immune response so variable?

-  What is the role of cytokine signaling in explaining immune response
   variability?

Some factors to be included in future developments of the model are:

-  How is the virus transported in the extracellular environment and
   mucus?

-  What is the role of the humoral immune response in controlling the
   viral spreading?

The model includes a representation of extracellular virus in the mucus,
epithelial cells and immune cells. It also includes the processes of
epithelial cell infection by extracellular virus, viral replication and
cell damage inside epithelial cells, release of viruses by epithelial
cells, immune cell response to infected epithelial cells and immune cell
killing of infected and non-infected epithelial cells. 

At the epithelial cell level, the model accounts for internalization,
replication, release and clearance of viral particles, as well as for
induced cell apoptosis by either viral damage or immune cytotoxicity.

-  **Viral internalization**: model of viral binding to cell receptors,
   endocytosis and release of genetic material into the cytoplasm 

-  **Viral replication**: model of replication of viral genetic
   material, transcription, translation and virion packing

-  **Viral release**: model of the release of newly assembled virions
   into the extracellular environment

-  **Viral damage:** model of accumulated damage to the cell due to
   viral load

-  **Cell death:** model of cell death due to accumulated damage from
   viral infection or by cytotoxicity from immune response

At the immune cell level, the model accounts for recruitment and
chemotaxis of immune cells due to cytokine signaling, the cytotoxic
effect on infected epithelial cells as well as the clearance of immune
cells.

-  **Immune cell recruitment**: model of immune cell recruitment and
   infiltration into the tissue by signaling molecules produced in
   response to viral insult on infected cells

-  **Immune cell chemotaxis**: model of immune cell movement guided by
   the difference in concentration of a signal represented as a
   chemical field

-  **Immune cell cytotoxicity**: model of cell-dependent cytotoxic
   effect of immune cells on infected cells

-  **Immune cell clearance**: model of immune cell-accumulated damage,
   cell death and clearance from the tissue

At the tissue level, the model accounts for the extracellular transport
of viral particles and can be extended to incorporate cytokine transport
and recovery by reepithelialization. 

-  **Viral transport**: model of diffusion and spreading of viral
   particles in the extracellular environment

-  **Cytokine transport**: model of transport of small immune signaling
   molecules in the extracellular environment (to be included)

-  **Tissue recovery**: model of recovery of the tissue by
   reepithelialization following cell death (to be included)

.. _fig1:

.. figure:: https://raw.githubusercontent.com/covid-tissue-models/covid-tissue-response-models/master/CC3D/Models/BiocIU/SARSCoV2MultiscaleVTM/media/image1.png
   :width: 4in
   :height: 3.39167in
   :align: center

   Conceptual Model

.. _model_description_end:

.. _model_implementation_start:

Model Implementation
====================

Epithelial cells can adopt one of three different phenotypes:
uninfected, infected and dead. Uninfected cells can absorb viral
particles from the extracellular environment but do not release newly
assembled particles until a critical viral load is reached. Once the
critical viral load is reached, uninfected cells change their cell type
to become infected cells. Infected cells both absorb and release viral
particles from the extracellular environment. Infected cells can become
uninfected cells by clearing their viral load. Infected cells can also
trigger apoptosis and become dead cells either by reaching a critical
viral load or by cytotoxic interaction with immune cells.

Intracellular state models
--------------------------

The submodels governing the intracellular state of uninfected and
infected cells associated with viral infection are implemented as
follows. 

Viral internalization
~~~~~~~~~~~~~~~~~~~~~

Uninfected and infected cells have the ability of absorbing diffusive
viral particles from the extracellular viral field. The uptake
probability
:math:`P\left( \text{Uptake}\left( \text{cell} \right) > 0 \right)` for each
cell occurs according to a Hill equation of the total amount of
diffusive viral particles in the domain of the cell
:math:`V\left( \text{cell} \right)`. 

.. math:: P\left( \text{Uptake} > 0 \right) = \frac{V\left( \text{cell} \right)^{h}}{V\left( \text{cell} \right)^{h} + V_{\text{half}}}

Where :math:`h` is a Hill coefficient and :math:`V_{\text{half}}` is the
local amount of the viral field at which
:math:`P\left( \text{Uptake} > 0 \right) = 0.50`. When uptake occurs, the
uptake rate is proportional to the local amount of the viral field and a
prescribed uptake efficiency :math:`r_{e}`, and saturates at a
prescribed threshold :math:`U_{\text{threshold}}`, 

.. math::

   \text{Uptake}  = \left\{ \begin{matrix}
   r_{e}V\left( \text{cell} \right)      & V\left( \text{cell} \right) < U_{\text{threshold}} \\
   U_{\text{threshold}} & V\left( \text{cell} \right) > U_{\text{threshold}} \\
   \end{matrix} \right.

The amount absorbed by each cell is subtracted from the viral field and
passed to the cell’s instance of the viral replication model according
to conservation of species.

Viral replication
~~~~~~~~~~~~~~~~~

A system of ordinary differential equations modeling the viral
replication process is assigned to each uninfected and infected cell.
The model contains four variables representing different states of the
viral replication process: unpacking :math:`U`, replicating :math:`R`,
packing :math:`P`,  and assembly of new virion capsids :math:`A`. 

.. math:: \frac{\text{d}U}{\text{d}t} = \text{Uptake} - r_{u}U

.. math:: \frac{\text{d}R}{\text{d}t} = r_{u}U + r_{\text{max}}\frac{R}{R + r_{\text{half}}} - r_{t}R

.. math:: \frac{\text{d}P}{\text{d}t} = r_{t}U - r_{p}P

.. math:: \frac{\text{d}A}{\text{d}t} = r_{p}P - \text{Secretion}

Here :math:`r_{u}` is the unpacking rate, :math:`r_{\text{max}}` is the
maximum replication rate, :math:`r_{t}` is the translation rate and
:math:`r_{p}` is the packing rate. The regulation of replication is
represented by a Michaelis-Menten function of the amount of replicating
viral material :math:`\frac{R}{R\  + r_{\text{half}}}`, where
:math:`r_{\text{half}}` is the amount of :math:`R` at which the
replicating rate is :math:`\frac{r_{\text{max}}}{2}`. The viral replication
model is specified as a readily sharable Antimony string that can be
implemented as a standalone using the Tellurium package. The number of
newly assembled virion capsids is passed to the cell’s instance of the
viral release model. 

Viral release
~~~~~~~~~~~~~

Infected cells have the ability to secrete diffusive viral particles
into the extracellular viral field. The total amount released is
proportional to the state variable for assembled virions from the viral
replication model. 

.. math:: \text{Secretion} = r_{s}A

Here :math:`r_{s}` is the secretion rate of viral particles. The amount
released by each cell is subtracted from the cell’s state variable for
assembled virions and passed to the source term of the extracellular
viral field according to conservation of species. 

Virally induced apoptosis
~~~~~~~~~~~~~~~~~~~~~~~~~

Each infected cell is assigned a survival probability. Once the state
variable for assembled virions from the viral replication model reaches
a prescribed critical threshold in a cell, the probability of cell
survival is evaluated against a uniformly distributed random variable.
Surviving cells remain infected and their survival is not re-evaluated.
Dying cells change cell type to dead cell and their instances of the
viral internalization, replication and release models are disabled. 

Immune response models
----------------------

Immune cells infiltrate the tissue and move up the gradient of the
extracellular viral field. The viral field is used as a proxy for
cytokines. Immune cells can induce cytotoxicity in infected cells and
trigger apoptosis. Immune cells are cleared out from the tissue.
Submodels of immune response are implemented as follows. 

Immune cell recruitment
~~~~~~~~~~~~~~~~~~~~~~~

Immune cells are seeded into the simulation space at a rate determined
by a prescribed seeding probability. At each simulation step the seeding
probability is evaluated against a uniformly distributed random
variable. To determine the seeding location, the simulation space is
randomly sampled, and immune cells are seeded at the unoccupied location
with the highest amount of the viral field. If no location is
unoccupied, then the immune cell is not seeded. 

Immune cell chemotaxis
~~~~~~~~~~~~~~~~~~~~~~

Immune cells experience a motile force as a response to a signaling
field. Currently, the viral field is used as a proxy of cytokine
signaling molecules. The chemotactic function measures the local
gradient of the viral field and computes the effective energy
:math:`E_{\text{chemotaxis}}` associated with the gradient according to
a prescribed chemotactic sensitivity parameter chemotaxis. The
chemotactic effective energy term is saturated by normalizing the
chemotactic sensitivity parameter by the local concentration
:math:`V\left( \text{cell} \right)`.

.. math:: E_{\text{chemotaxis}} = \frac{\lambda_{\text{chemotaxis}}\nabla V}{1 - V\left( \text{cell} \right)}

Immune cell cytotoxicity
~~~~~~~~~~~~~~~~~~~~~~~~

Immune cells kill infected cells by direct contact. At each simulation
step, neighbors of infected cells are evaluated. Apoptosis is triggered
in an infected cell if it has an immune cell as one of its neighbors.
The infected cell changes its cell type to dead cell and its instances
of the viral internalization, replication and release models are
disabled.

Immune cell clearance
~~~~~~~~~~~~~~~~~~~~~

Each infected immune cell is assigned a dying probability. For each
simulation step, the dying probability is evaluated against a uniformly
distributed random variable for every infected cell. Clearance is
achieved by setting the immune cell volume constraint to zero.

Transport models
----------------

The extracellular viral field is used to represent the transport of
viral particles across the tissue over time. Rates of secretion into the
viral field are determined by the output of the viral release model.
Rates of absorption from the viral field are determined by the viral
internalization model.

Viral transport
~~~~~~~~~~~~~~~

The change in concentration of the viral field at each location is
calculated using a partial differential equation solver of a
reaction-diffusion equation. 

.. math:: \frac{\partial V\left( x \right)}{\partial t} = D\mathrm{\Delta}V - cV\left( x \right) - \text{Uptake}\left( \text{Cell}\left( x \right) \right) + \text{Secretion}\left( \text{Cell}\left( x \right) \right)

Transport parameters such as the diffusion constant :math:`D` and decay
rate :math:`c` are estimated from the literature. Conversion factors are
used to translate experimental parameter values to internal simulation
parameters.

.. _fig2:

.. figure:: https://raw.githubusercontent.com/covid-tissue-models/covid-tissue-response-models/master/CC3D/Models/BiocIU/SARSCoV2MultiscaleVTM/media/image2.png
   :width: 5in
   :height: 2.52014in
   :align: center

   Interactions in the Tissue Model

.. _fig3:

.. figure:: https://raw.githubusercontent.com/covid-tissue-models/covid-tissue-response-models/master/CC3D/Models/BiocIU/SARSCoV2MultiscaleVTM/media/image3.png
   :width: 5in
   :height: 1.80833in
   :align: center

   Interactions in the Viral Replication Model

.. _model_implementation_end: