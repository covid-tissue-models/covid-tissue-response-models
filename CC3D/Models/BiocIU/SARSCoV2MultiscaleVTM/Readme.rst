
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
