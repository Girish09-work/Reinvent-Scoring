+-----------------------------------+-----------------------------------+
|                                   | *Briefings in Bioinformatics*,    |
|                                   | 2023, **24(5)**, 1--13            |
+===================================+===================================+
|                                   | > **https                         |
|                                   | ://doi.org/10.1093/bib/bbad323**\ |
|                                   | > Advance access publication date |
|                                   | > 5 September 2023 **Problem      |
|                                   | > Solving Protocol**              |
+-----------------------------------+-----------------------------------+

> **3D based generative PROTAC linker design with reinforcement
> learning**

+-----------------------+-----------------------+-----------------------+
| Baiqing Li            | ![](vert              | > , Ting Ran and      |
|                       | opal_6f7e04f02c2f4d50 | > Hongming Chen       |
|                       | 95038d2511b10a73/medi |                       |
|                       | a/image1.png){width=" |                       |
|                       | 0.1527777777777778in" |                       |
|                       | height="0             |                       |
|                       | .1527777777777778in"} |                       |
+=======================+=======================+=======================+
+-----------------------+-----------------------+-----------------------+

> Corresponding author. Hongming Chen, Guangzhou National Laboratory,
> Guangzhou 510005, Guangdong Province, China. Tel.: +86-19949451768;
> Fax: +86 020-84098081; E-mail: chen_hongming@gzlab.ac.cn
>
> Abstract
>
> Proteolysis targeting chimera (PROTAC), has emerged as an effective
> modality to selectively degrade disease-related proteins by harnessing
> the ubiquitin-proteasome system. Due to PROTACs' hetero-bifunctional
> characteristics, in which a linker joins a warhead binding to a
> protein of interest (POI), conferring specificity and a E3-ligand
> binding to an E3 ubiquitin ligase, this could trigger the
> ubiquitination and transportation of POI to the proteasome, followed
> by degradation. The rational PROTAC linker design is challenging due
> to its relatively large molecular weight and the complexity of
> maintaining the binding mode of warhead and E3-ligand in the binding
> pockets of counterpart. Conventional linker generation method can only
> generate linkers in either 1D SMILES or 2D graph, without taking into
> account the information of ternary structures. Here we propose a novel
> 3D linker generative model PROTAC-INVENT which can not only generate
> SMILES of PROTAC but also its 3D putative binding conformation coupled
> with the target protein and the E3 ligase. The model is trained
> jointly with the RL approach to bias the generation of PROTAC
> structures toward pre-defined 2D and 3D based properties. Examples
> were provided to demonstrate the utility of the model for generating
> reasonable 3D conformation of PROTACs. On the other hand, our results
> show that the associated workflow for 3D PROTAC conformation
> generation can also be used as an efficient docking protocol for
> PROTACs.
>
> Keywords: PROTAC linker design, generative model, reinforcement
> learning
>
> INTRODUCTION
>
> Protein hydrolysis targeting chimeric (PROTAC) is emerging as a
> promising technology in the realm of drug discovery. Unlike
> traditional occupancy-driven drug design, it aims to degrade drug
> target protein through hijacking the ubiquitin proteasome system (UPS)
> in eukaryotic cells \[1, 2\]. A standard PROTAC molecule consists of
> three parts, namely, a warhead fragment that binds to the protein of
> interest (POI), a warhead fragment that binds
>
> may lead to drug resistance caused by mutations at the binding site of
> POI \[4, 5\]. Since the concept of PROTAC was proposed by Crews *et
> al.* \[1, 7\], the great potential of the technology has attracted
> huge interests among researchers from academics as well as
> pharmaceutical industry and similar ideas has also been applied in
> other areas \[8--13\]. Currently, a few PROTAC drugs have moved into
> clinical stage \[14\]. Some recent reviews have comprehensively
> summarized the characteristics of PROTACs
>
> to the E3 ubiquitin ligase (E3 ligand) and a linker segment that \[4,
> 5\].
>
> connects the two fragments. This assembly method allows the PROTAC
> molecule to bring ubiquitinase to the vicinity of POI to ubiquitinate
> the POI and the ubiquitinated POI is then degraded by the UPS. It has
> exhibited some unique advantages over traditional small-molecule
> drugs. For example, as UPS is involved in degra-dation of *\>*80% of
> proteins in cells \[3\], the PROTAC technology is theoretically able
> to degrade most proteins in cells. Second, PROTAC mediated degradation
> process does not rely on high affinity of the warhead to POI \[4\], so
> there is less restriction on the binding site of POI. This greatly
> expands drug target space so that many traditionally undruggable
> proteins may become druggable via PROTAC modalities \[5, 6\]. In
> addition, PROTAC molecule can eliminate all functions of POI through
> protein degradation, rather than just blocking warhead mediated
> biological functions, which
>
> So far, thousands of PROTAC molecules have been reported and the
> PROTAC-DB database \[15\] curated *\>*2000 published PROTAC molecules.
> Compared with the chemical space of traditional drug-like small
> molecule, the PROTAC space is largely different due to their
> relatively large molecular size and still under-represented. So far,
> up to 60 E3 ligases have been found, but existing E3 ligands of PROTAC
> are mainly focused on the ligand set of von Hippel--Lindau and
> cereblon (CRBN) ligases \[16\]. The degrada-tion ability of most E3
> ligases has not been fully exploited. On the other hand, studies show
> that the linker segment is also critical for PROTAC design. It has
> been shown that the linker segment is closely associated with the
> pharmacokinetic proper-ties and protein degradation ability of PROTACs
> \[16, 17\]. Espe-cially, it plays a pivotal role in maintaining the
> conformational
>
> **Baiqing Li** is currently a Research Associate in Guangzhou National
> Laboratory. His research interests lie in the deep learning, knowlege
> graph embedding, drug discovery and computational biology.
>
> **Ting Ran** is currently an Associate Professor in Guangzhou National
> Laboratory. His research interests lie in the drug discovery,
> computer-aided drug design and machine learning.
>
> **Hongming Chen** has more than 20 years experience of pharmaceutical
> industry, working as computational chemist. He currently holds a
> Principal Investigator position of Guangzhou National Laboratory at
> Guangzhou. His research group concentrates on method development of
> computational chemistry, especially on AI augmented drug design.
>
> **Received:** June 5, 2023. **Revised:** August 6, 2023. **Accepted:**
> August 20, 2023\
> © The Author(s) 2023. Published by Oxford University Press. All rights
> reserved. For Permissions, please email: journals.permissions@oup.com

+-----------------+-----------------+-----------------+-----------------+
| > 2             | \|              | > *Li* et al.   |                 |
+=================+=================+=================+=================+
| > stability of  |                 |                 | > METHODS       |
| > PROTAC        |                 |                 |                 |
| > ternary       |                 |                 |                 |
| > structure     |                 |                 |                 |
| > (PTS) (ie.    |                 |                 |                 |
| > the           |                 |                 |                 |
| > POI-PROTAC-   |                 |                 |                 |
+-----------------+-----------------+-----------------+-----------------+

> E3 complex), since an appropriate linker fragment can ensure the
> energetically stable binding conformation of E3 ligand and warhead of
> POI \[17\], thus forming a stable ternary complex. However, available
> experience for linker design is still limited. Most of earlier linker
> fragments come from derived PEG chain due to its high flexibility and
> ease of synthesis. On this basis, various flexible linker fragments
> have been developed for PRO-TAC. Recently, relatively rigid linker
> fragments have also been proposed exemplified by the androgen receptor
> targeted PROTAC compounds in the clinical stage \[18\]. Thus,
> developing method-ology for PROTAC linker design is urgently needed
> for PROTAC discovery.
>
> Recently, deep learning based molecular generative model has shown
> significant advancement \[19--24\] and have been applied on fragment
> linking which could be used for fragment-based drug design and PROTAC
> design. For example, SyntaLinker \[25\], a SMILES based linker
> generative model based on the well-known transformer architecture, was
> recently been proposed to design linker to fuse two terminal fragments
> into a complete molecule. Zheng *et al.* \[26\] recently has applied a
> similar methodology on PROTAC design. To overcome the deficiency of
> PROTAC molecules for model training, a large collection of
> quasi-PROTAC molecules that have similar characteristics to PROTACs
> was used to pre-train a prior model, which was then fine-tuned by
> using actual PRO-TACs augmented with randomized SMILES of fragments.
> More-over, they utilized reinforcement learning (RL) to optimize the
> model to generate PROTACs with potentially better pharmacoki-netic
> properties. This model showed superior performance on PROTAC
> generation over other linker generative models. In con-trast, Delinker
> \[27\] is a graph-based 2D linker generative model and it was claimed
> that it could be used to design PROTAC molecules. Igashov *et al.*
> \[28\] proposed Difflinker model by inte-grating diffusion model with
> information of protein pocket for linker generation, and multiple
> fragments can be assembled into a complete molecule. Guo *et al.*
> \[29\] reported Link-INVENT, a 2D linker generative model which
> demonstrated potential to gener-ate both small molecules that conform
> to Lipinski 'rule of five' \[30\] and PROTACs. The algorithm was
> essentially an expansion of the widely used REINVENT model \[31\]
> incorporating a reward func-tion to optimize the length, linearity and
> flexibility of generated linkers.
>
> Although these progresses in the development of PROTAC gen-erative
> model, available methods are mainly focused on the gen-eration of 2D
> PROTAC-like structures and did not consider the fea-sibility of a
> PROTAC fitting in the binding site of the PTS complex. To address this
> concern, we proposed a novel linker generative model, named as
> PROTAC-INVENT, utilizing three-dimensional PTS data as constraints.
> The workflow consists of two stages: generation of PROTAC linker
> structures (in SMILES format) from normal generative model and then
> creation of 3D binding con-formation of PROTAC under the constraint of
> the binding pocket of a reference PTS. The model was fine-tuned via RL
> and gener-ated PROTACs were scored based on their docking poses in the
> PTS binding site. Compared with other linker generative models, one
> unique feature of PRO-INVENT is that a practical workflow was
> developed to enable generation of 3D binding conformation of PROTAC in
> the binding site with relatively low computation cost, so that
> generation of PROTAC molecules is optimized by considering its fitness
> with the PTS binding site. Thus, this method is expected to be more
> effective in the generation of PROTAC molecules which need to satisfy
> the 3D requirement of PTS bind-
>
> **Model overview**
>
> PROTAC-INVENT takes a pair of fragment (E3-ligand and warhead for POI)
> and a reference PTS as input, and returns generated link-ers and the
> predicted binding conformations of the full PROTAC compounds in the
> binding site of reference PTS (as shown in Figure 1). T wo modules are
> integrated in PROTAC-INVENT. First, a generative model is trained to
> produce chemically feasible linker SMILES strings, following by
> forming complete PROTAC molecules together with the specified warhead
> and E3-ligand. REINVENT was used as the core generative engine for
> linker design and it was previously developed as an efficient
> generative model for *de novo* structure design. RL was utilized to
> search chemical space for optimization of molecular properties.
> Second, a 3D structure-based workflow was implemented to generate
> docking conforma-tions of PROTACs and provide scores accordingly. The
> 3D based scores were then utilized to drive the RL (as shown in Figure
> 2).
>
> **Generation of PROTAC 3D conformations**
>
> Previously, Link-Invent was designed to generate PROTAC link-ers via
> RL, but their scoring function was based on purely 2D related metrics.
> Here we implemented a practical workflow to enable predicting 3D
> binding conformations of generated PROTAC molecules and score them by
> taking into account the structure information of a reference PTS. It
> takes four steps for the protocol to convert the SMILES of PROTAC to
> 3D binding conformation given a reference PTS as shown in Figure 3:
>
> \(1\) Generate initial conformers of PROTAC in the vicinity of
> reference ligand
>
> First of all, the SMILES of linker was generated by REINVENT prior
> model and then merged with the supplied warhead pair (as part of the
> input) to form the SMILES of full PROTAC. *Omega* \[32\] program was
> then utilized to convert SMILES to a set of initial 3D conformations.
> The 3D conformations were then super-imposed with the bound
> conformation of the reference PROTAC by using commercially available
> *ROCS* program \[33\]. The Com-boScore which is a weighted sum of
> shape and pharmacophore similarity was used to measure similarity
> between the reference conformer and the aligned conformation. Due to
> the large size of PROTAC, this superimposing is usually not able to
> create good alignment but still can serve the purpose of bringing the
> linker part of the PROTAC to the vicinity of the reference linker.
>
> \(2\) Merge the linker of PROTAC with the warheads of reference ligand
>
> Once the alignment was done, the warheads of the generated PROTAC were
> removed and only the linker atoms were remained. The warheads of
> reference ligand were then copied and merged with the linker fragment
> of the generated PROTAC to form a re-combined PROTAC conformation
> (RPC). After this operation, the coordinates of warhead atoms in the
> RPC were exactly the same to that of reference ligand, but the bond
> lengths, dihedral angles between terminal atoms of the linker and
> warheads were not corrected yet and need to be fixed.
>
> \(3\) Optimize the conformation of RPC
>
> In order to optimize the conformation of RPC while keeping its warhead
> parts as close to that of the reference PROTAC con-formation as
> possible, a constrained molecular mechanics opti-mization was then
> carried out. Here we employed the MacroModel
>
> ing site. module of*Schrödinger* software \[34\] to achieve the
> constrained

  -----------------------------------------------------------------------
  *3D based generative    \|                      3
  PROTAC linker design*                           
  ----------------------- ----------------------- -----------------------

  -----------------------------------------------------------------------

![](vertopal_6f7e04f02c2f4d5095038d2511b10a73/media/image2.png){width="5.0in"
height="2.279166666666667in"}

> Figure 1. Illustration of PROTAC-INVENT model. The model is intended
> to generate PROTACs whose E3-ligand and POI warhead motifs can adopt
> similar binding conformation as that of the PROTAC in the reference
> PTS.

![](vertopal_6f7e04f02c2f4d5095038d2511b10a73/media/image3.png){width="5.001388888888889in"
height="3.9152766841644793in"}

> Figure 2. The general workflow of PROTAC-INVENT.
>
> minimization by imposing constraints on warhead atoms of RPC. The
> specific implementation, parameter setting and precautions can be
> found in supporting materials. After 200 iterations, the minimized RPC
> was saved as the optimized PROTAC conformation
>
> due to the large search space, therefore it is not guaranteed that, in
> docking pose, the POI warhead and E3-ligand motifs can adopt similar
> poses as in the reference ligand. To address this problem, we choose
> to use the 'local-only' mode (or refine
>
> (OPC). mode) to do docking, in which the initial pose sampling stage
>
> \(4\) Docking OPC into the binding pocket of PT[S]{.underline}
>
> To further evaluate generated PROTAC, its OPCs were docked into the
> binding pocket of PTS, which is composed by E3 ligase and POI, to
> obtain their docking scores. The docking procedure usually contains
> two parts: initial conformation/pose sampling to obtain multiple
> starting points and the following conforma-tion minimization from the
> sampled starting points. For PRO-TAC molecule, its large size will
> result in large conformational search space and accordingly lead to
> long computation time.
>
> was skipped and the input conformation was taken as the sole starting
> point for optimization. The 'local-only' docking solution was taken as
> the final conformation of PROTAC molecule and its docking score was
> used as a scoring component for RL. The docking component from
> DockStream \[35\] that is fully compatible with REINVENT \[31\] was
> served as the interface between REIN-VENT and docking software and
> DockStream supports various docking programs. In this work, we chose
> *AutoDock Vina* \[36--38\] for doing 'local-only' docking; its docking
> score was used as a scoring component to evaluate the PROTAC molecule
> in RL
>
> Searching for global minimum of PROTAC also becomes difficult loops.

+-----------------------+-----------------------+-----------------------+
| > 4                   | \|                    | > *Li* et al.         |
+=======================+=======================+=======================+
+-----------------------+-----------------------+-----------------------+

![](vertopal_6f7e04f02c2f4d5095038d2511b10a73/media/image4.png){width="6.5in"
height="6.5111100174978125in"}

> Figure 3. The workflow of generating and scoring 3D binding
> conformation of PROTAC.

+---------+---------+---------+---------+---------+---------+---------+
| > **    | >       | *S(x)*  | *i*     | *iwi* × | > *iwi* | \(1\)   |
| Scoring | number. | =       |         | *pi(x)  |         |         |
| > of    |         |         |         | iwi*    |         |         |
| > PR    |         |         |         |         |         |         |
| OTACs** |         |         |         |         |         |         |
+=========+=========+=========+=========+=========+=========+=========+
| > Mo    |         |         |         |         |         |         |
| lecular |         |         |         |         |         |         |
| > gen   |         |         |         |         |         |         |
| erative |         |         |         |         |         |         |
| >       |         |         |         |         |         |         |
|  models |         |         |         |         |         |         |
| > such  |         |         |         |         |         |         |
| > as    |         |         |         |         |         |         |
| > R     |         |         |         |         |         |         |
| EINVENT |         |         |         |         |         |         |
| >       |         |         |         |         |         |         |
| utilize |         |         |         |         |         |         |
| > RL to |         |         |         |         |         |         |
+---------+---------+---------+---------+---------+---------+---------+
| > guide |         |         |         |         |         |         |
| > st    |         |         |         |         |         |         |
| ructure |         |         |         |         |         |         |
| > gene  |         |         |         |         |         |         |
| ration. |         |         |         |         |         |         |
| > In    |         |         |         |         |         |         |
| > RL,   |         |         |         |         |         |         |
| > co    |         |         |         |         |         |         |
| mpounds |         |         |         |         |         |         |
| > are   |         |         |         |         |         |         |
| > first |         |         |         |         |         |         |
| >       |         |         |         |         |         |         |
|  scored |         |         |         |         |         |         |
+---------+---------+---------+---------+---------+---------+---------+
|         |         |         |         | 1*/*    |         |         |
+---------+---------+---------+---------+---------+---------+---------+
| > and   |         | *S(x)*  |         |         |         | \(2\)   |
| > the   |         | =       |         |         |         |         |
| >       |         |         |         |         |         |         |
|  scores |         |         |         |         |         |         |
| > are   |         |         |         |         |         |         |
| > then  |         |         |         |         |         |         |
| > fed   |         |         |         |         |         |         |
| > back  |         |         |         |         |         |         |
| > to    |         |         |         |         |         |         |
| >       |         |         |         |         |         |         |
|  neural |         |         |         |         |         |         |
| > n     |         |         |         |         |         |         |
| etworks |         |         |         |         |         |         |
| > to    |         |         |         |         |         |         |
| >       |         |         |         |         |         |         |
|  update |         |         |         |         |         |         |
+---------+---------+---------+---------+---------+---------+---------+
|         |         |         |         | > *p    |         |         |
|         |         |         |         | i(x)wi* |         |         |
+---------+---------+---------+---------+---------+---------+---------+
| > their |         |         |         |         |         |         |
| > par   |         |         |         |         |         |         |
| ameters |         |         |         |         |         |         |
| > to    |         |         |         |         |         |         |
| > i     |         |         |         |         |         |         |
| ncrease |         |         |         |         |         |         |
| > the   |         |         |         |         |         |         |
| > prob  |         |         |         |         |         |         |
| ability |         |         |         |         |         |         |
| > of    |         |         |         |         |         |         |
| > str   |         |         |         |         |         |         |
| uctures |         |         |         |         |         |         |
| > with  |         |         |         |         |         |         |
+---------+---------+---------+---------+---------+---------+---------+
| >       |         |         |         |         |         |         |
|  better |         |         |         |         |         |         |
| >       |         |         |         |         |         |         |
|  score. |         |         |         |         |         |         |
| > In    |         |         |         |         |         |         |
| > RE    |         |         |         |         |         |         |
| INVENT, |         |         |         |         |         |         |
| > two   |         |         |         |         |         |         |
| > types |         |         |         |         |         |         |
| > of    |         |         |         |         |         |         |
| > m     |         |         |         |         |         |         |
| ulti-co |         |         |         |         |         |         |
| mponent |         |         |         |         |         |         |
| >       |         |         |         |         |         |         |
| scoring |         |         |         |         |         |         |
+---------+---------+---------+---------+---------+---------+---------+

> function were defined, in which individual components of the scoring
> function can be either combined as a weighted sum or as a weighted
> product as shown in Equations 1 and 2. Given a sequence x, the weight
> coefficient of individual score component *pi(x)* reflects its
> importance in the overall score *S*(*x*). For the weighted sum score
> (Equation 1), the weight for each component should be a floating
> number in the \[0, 1\] range and the sum of all
>
> *Component of docking score*
>
> In current study, several scoring components were used for
> con-structing the over-all score of PROTACs. One component *Pd* is the
> docking score of PROTAC and it was transformed to a floating number
> the \[0, 1\] range according to Equation 3, which is a reversed
> sigmoid function \[31\]:

+---------+---------+---------+---------+---------+---------+---------+
| >       | *P(x)*  | > 1 10  | *kv     |         | 10      | \(3\)   |
| weights | 1*/*    |         | high* + |         |         |         |
| >       |         |         | *low*   |         |         |         |
|  should |         |         |         |         |         |         |
| > be    |         |         |         |         |         |         |
| > 1.0.  |         |         |         |         |         |         |
| > For   |         |         |         |         |         |         |
| > the   |         |         |         |         |         |         |
| > w     |         |         |         |         |         |         |
| eighted |         |         |         |         |         |         |
| >       |         |         |         |         |         |         |
| product |         |         |         |         |         |         |
| >       |         |         |         |         |         |         |
|  score, |         |         |         |         |         |         |
| > the   |         |         |         |         |         |         |
| >       |         |         |         |         |         |         |
| overall |         |         |         |         |         |         |
+=========+=========+=========+=========+=========+=========+=========+
|         |         |         |         |         |         |         |
+---------+---------+---------+---------+---------+---------+---------+
| > score | *d* =   | > \+ ∗∗ | > ∗ −   | ∗       | *high*  |         |
| > is    |         |         |         |         | −*low*  |         |
| >       |         |         |         |         |         |         |
| defined |         |         |         |         |         |         |
| > as    |         |         |         |         |         |         |
| > E     |         |         |         |         |         |         |
| quation |         |         |         |         |         |         |
| > 2 and |         |         |         |         |         |         |
| > the   |         |         |         |         |         |         |
| >       |         |         |         |         |         |         |
| weights |         |         |         |         |         |         |
| > can   |         |         |         |         |         |         |
| > be    |         |         |         |         |         |         |
| > any   |         |         |         |         |         |         |
| > f     |         |         |         |         |         |         |
| loating |         |         |         |         |         |         |
+---------+---------+---------+---------+---------+---------+---------+

  -----------------------------------------------------------------------
  *3D based generative    \|                      5
  PROTAC linker design*                           
  ----------------------- ----------------------- -----------------------

  -----------------------------------------------------------------------

![](vertopal_6f7e04f02c2f4d5095038d2511b10a73/media/image5.png){width="6.501388888888889in"
height="2.9833333333333334in"}

> Figure 4. The Ps and raw ADV score of two generated PROTAC examples.
> The green structures are docking poses and the gray structures are
> reference poses.
>
> Here *v* is the docking score (more negative is better). *k* value is
> an empirical value that represents the steepness of score function
> curve (the default value is 0.25 here). Parameters *low, high* refer
> to the low and high end cut-off values of docking score and can be
> adjusted by user. Naturally, their exact values depend on different
> ternary complex systems. Taking the BTK PTS system for example, when
> the reference PROTAC molecule was re-docked into PTS binding site
> using Autodock-Vina, its docking score was −10.9, *low, high* values
> were set as −15 and −5, representing most favorable and unfavorable
> docking scores respectively.
>
> length which is the shortest 2D bond distance between warheads. If the
> linker length of generated structure is within the range, the score is
> 1.0. Otherwise, the score is 0.0.
>
> *Component for the number of aromatic ring*
>
> Studies have shown \[39, 40\] that multiple aromatic rings in linker
> of PROTAC molecules could be detrimental to the PK proper-ties of
> PROTACs. PROTAC-INVENT can customize the number of aromatic ring in
> the generated linker. Here, the default allowed number of aromatic
> ring is 1. If the number of aromatic ring of linker is less than or
> equal to the cut-off value, this score is 1.
>
> *Component of post-dock shape similarity* Otherwise, the score is 0.
>
> Another component, *Ps*, with a value in the range of \[0,1\], was set
> to measure the shape similarity between warheads (including the
> ligands for E3 ligase and POI) of docking conformation and that of the
> reference ligand. After 'local only' docking of DPC in the binding
> site of PTS, sometimes there is a large positional shift for the
> warheads of PROTACs comparing to the warheads of reference PROTAC. The
> larger the *Ps* value is, the more similar the warheads of PROTAC are
> to the reference. Here,*ROCS* was used to only calculate the shape
> similarity on site and no alignment optimization was done. For example
> in Figure 4, two PROTAC molecules *a* and *b* were generated by prior
> model. After docking procedure, the *Ps* scores of molecule *a* and
> *b* were 0.903 and 0.687, respectively, while their docking scores
> were roughly equal (ie. -11.37 versus −11.29). Clearly, in the docking
> poses, the warheads of molecule *a* had less deviation to those of
> reference ligand than molecule *b*.
>
> *Component for substructure alert*
>
> To clean the generated structures from generative model, a set of 24
> unwanted substructures (as shown in in supporting information) was
> defined. If the generated linker contains these substructures, the
> alert score is 0, otherwise the score is 1. The alert substructures
> can be user specified.
>
> *Component for linker length*
>
> The linker length is critical for maintaining proper conformation of
> PROTACs in PTS. PROTAC-INVENT can customize the range of linker length
> of PROTAC through a score component of linker
>
> **Computational details and datasets**
>
> REINVENT model was used as the core engine of PROTAC-INVENT to
> generate linkers in SMILES format. The relevant source code was mainly
> adapted from Link-INVENT \[29\], an extension of REINVENT model for
> linker design by Guo *et al*. Details should be referred to the
> original literature. In this work, the prior model of Link-INVENT was
> used as our prior model. This model was originally trained on a
> dataset processed from ChEMBL database, which comprises 10 313 109
> structure triplets which are formed as ('Warhead1 \| Warhead2, Linker,
> Full Molecule'). The processing details should be referred to the
> original paper \[29\]. During RL, the training epoch was set to 200
> and batch size was set to 100. Other parameters and high-parameters of
> training process can be found in supporting materials. In current
> study, Bruton's tyrosine kinase (BTK) PROTAC was chosen as the test
> case. Published PTS crystal structure 6W8I \[41\] was chosen as the
> reference structure for structure generation. For the diversity filter
> setting \[29, 31, 42, 43\] in the scoring function of RL, the
> 'IdenticalMurckoScaffold'option was used to increase the diversity of
> generated linkers. The'minscore' and 'minsimilarity' paremeters were
> all set as 0.0, so that all generated PROTACs during RL run will be
> kept.
>
> RESULT AND DISCUSSION\
> **Design of BTK PROTAC through PROTAC-INVENT model**
>
> To validate our methodology, design of Bruton's tyrosine kinase (BTK)
> PROTAC was selected for our case study. BTK is non-receptor

+-----------------------+-----------------------+-----------------------+
| > 6                   | \|                    | > *Li* et al.         |
+=======================+=======================+=======================+
+-----------------------+-----------------------+-----------------------+

![](vertopal_6f7e04f02c2f4d5095038d2511b10a73/media/image6.png){width="6.501388888888889in"
height="2.9499989063867016in"}

> Figure 5. The BTK-degrader-cIAP1 ternary complex crystal structure.
> The red box of PROTAC molecules represents the linker part and the two
> terminal fragments are POI warhead and E3-ligand, respectively.

![](vertopal_6f7e04f02c2f4d5095038d2511b10a73/media/image7.png){width="6.5in"
height="2.5722222222222224in"}

> Figure 6. The training curves based on two aggregation method of score
> components at different link length. The changes of average total
> score using (**a**) product aggregation and (**b**) additive
> aggregation.
>
> tyrosine kinase that plays a role in the maturation of B cells \[44\].
> It represents a validated drug target in leukemia and lymphoma
> \[45--47\] with well-described potent and selective ligands
> \[48--50\]. Within the protein degrader field, development of BTK
> PROTAC was reported, owing in part to observations related to its high
> apparent degradability \[51\]. Schiemer *et al.* \[41\] recently
> reported two BTK PTSs (PDB code: 6W8I/6W7O), demonstrating two
> distinct binding poses between the target and E3 ligase. In one PTS,
> the linker part is five repeated PEG units and the effective linker
> length is 15, and, for the other PTS, a pyrazine ring was used instead
> (the linker effective length is 7).
>
> We took one of the published BTK PTS (6W8I) as the refer-ence
> structure for linker generation. To explore the effect of the linker
> length on the generation of PROTAC molecules, we imposed constraint on
> linker length which refers to the topological bond distance between
> POI warhead and E3 ligand. Given the published BTK PTSs, we set
> several cut-offs for linker length, that is linker length cut-off
> (7--9), (7--11), (7--13) and (7--15), and the definition of linker
> score can be seen in previous section.
>
> **Analysis of different aggregation methods**\
> As described in previous section, two aggregation methods, pro-duction
> and summation form, were used for scoring molecules and their
> influence on training process were examined. Experi-ments with
> different aggregation method under different linker length constraints
> were carried out and their training curves are shown in Figure 6. When
> summation score was used, as shown in Figure 6b, during the training
> of 200 epochs, the curves of average total score were largely
> fluctuated in all models constrained with linker length, while the
> curves for production score (as shown in Figure 6a) behaved much
> better. In most of cases, the total score increased to around 0.8 in
> 50 epochs and can be stabilized more or less at the level for the rest
> of training. The results show that the training of production based
> aggregation method is easier than the summation based method. So in
> the following experiments, only production scoring models were
> discussed.
>
> As the total score contains multiple components, the weights of
> component might also influence the training. Among all the score
> components, ADV score (the docking score of *Autodock Vina*)

  -----------------------------------------------------------------------
  *3D based generative    \|                      7
  PROTAC linker design*                           
  ----------------------- ----------------------- -----------------------

  -----------------------------------------------------------------------

![](vertopal_6f7e04f02c2f4d5095038d2511b10a73/media/image8.png){width="6.501388888888889in"
height="2.5388888888888888in"}

> Figure 7. The training curves at different docking weight. The changes
> of (**a**) average total score and (**b**) ADV score.
>
> was regarded as the most important component for evaluating generated
> molecules in PTS. We paid more attention on this component and several
> weight values of ADV score were explored, while the weight values for
> other components were set as 1.0. The training curves are shown in
> Figure 7. It can be seen from Figure 7a that for all ADV weights, the
> average total score can quickly converge to 0.8. In the case of ADV
> weight of 4.0, the fluctuation of total score becomes larger than
> other weight level. The actual ADV scores during training can be seen
> in Figure 7b. It seems that the average docking score quickly reaches
> 0.7 and fluctuates between 0.7 and 0.9 for the rest of training
> epochs. There is no clear trend on how increasing weight of docking
> score would influence the docking score itself. This might due to the
> fact that slight change of the linker structure could change the
> docking score of the structure drastically. However, it can be seen
> that, in general, the total score can be improved during the RL.
>
> During RL, all generated PROTACs will be saved as poten-tial solutions
> for later inspection. We analyzed the number of molecules in the
> collection under different constraints of linker length and the result
> is shown in Figure 8. In this case, during the RL run, only constraint
> of link length was varied and all other parameters were the same. As
> shown in Figure 8, if no constraint of link length was imposed, the
> link length of generated solutions covers a very wide range (2--35),
> while adding constraint on linker length, the linker length of
> generated solutions is in the range of 5--15. As the range of
> constraint of link length becomes wider, the distribution of linker
> length of collected solutions also becomes wider, although it is not
> dramatically. This suggested that in RL framework, the restraining of
> link length does work.
>
> The effect of RL was also examined by looking at ratio of compounds
> which have good scores as shown in Table 1. As a comparison, without
> running RL, we randomly sampled the prior model (the model from
> Link-INVENT) to generate 20 K PROTACs as the baseline model and
> removed duplicates. For the ratio of compound whose total score is
> *\>* 0.8, the 'No RL' procedure is roughly 0.51, which is lower than
> that of other RL procedures. If we look into the ratio of *Vina*
> docking score (ie. ratio for ADV score *\<* −10.9, which is the
> docking score of reference PROTAC), a component of the total score,
> the difference is more significant. For PROTACs generated from the 'No
> RL' procedure, its ratio is about 0.35, while that of other RL
> procedures are *\>*0.54. Similar trend was observed for the ratio of
> *Ps\>* 0.9. This result shows
>
> ![](vertopal_6f7e04f02c2f4d5095038d2511b10a73/media/image9.png){width="3.0in"
> height="2.2694433508311462in"}
>
> Figure 8. Distribution of linker length of generated PROTACs in
> various RL runs.
>
> that the RL approach can generated more desirable solutions than
> sampling from prior model alone without running RL. Some example BTK
> PROTAC compounds and their potential docking conformations are
> illustrated in Figure 9. Additionally, the same protocol was applied
> on generating PROTACs for BAF \[52\] and BRD4 \[53\]. Their results
> and illustrative examples were included in the Supporting material.
>
> **The analysis of molecular dynamics simulations and MM-GBSA
> calculation**
>
> To further validate those generated PROTAC structures, molecular
> dynamic (MD) simulations and following molecular mechanics with
> generalized Born and surface area solvation (MM-GBSA) calculation were
> applied to estimate the binding affinities of some generated
> compounds. Taking BTK system (PDB code: 6W8I) as example, and manually
> selected 25 generated PROTACs, alone with the reference PROTAC
> structure in 6W8I, we performed MD simulations to capture the movement
> of mediated ternary complexes. Each MD simulation run for 50 ns and
> the last 5 ns (100 conformations) were extracted to calculate the
> MM-GBSA energy (after removing all water molecules). The 100 energy
> values were averaged to obtain the final binding affinity. The result
> is shown in in supporting information. It is shown that 5 of the 25
> selected PROTACs demonstrate better binding affinity than the

+-----------------------+-----------------------+-----------------------+
| > 8                   | \|                    | > *Li* et al.         |
+=======================+=======================+=======================+
+-----------------------+-----------------------+-----------------------+

![](vertopal_6f7e04f02c2f4d5095038d2511b10a73/media/image10.png){width="6.5in"
height="6.283333333333333in"}

Figure 9. Four generated PROTAC compounds for BTK system and their
docking conformations. PROTAC X-ray conformation was used as the
reference.

> reference PROTAC structure. Theoretically, these results provide some
> additional evidence on the reliability of our model. The MD study was
> conducted by using DESMOND \[54\] (Schrödinger release 2020-1). The
> GBSA calculations were performed using the Schrodinger Thermal MM-GBSA
> tool.
>
> **The performance of derived PROTAC docking protocol**
>
> Besides the generation of linker structures, the workflow for
> generation of PROTAC docking pose can also be used as a prac-tical
> docking protocol. To the best of our knowledge, there is still no
> docking software dedicated for PROTAC molecules. Given the large size
> of PROTAC compound, it would be a challenge to properly dock designed
> PROTAC compounds into the binding site of PTS. The PROTAC docking
> protocol used in PROTAC-INVENT can quickly generate docking pose of
> new PROTACs by mimicking the pose of reference PROTAC structure. To
> evaluate the docking performance of PROTAC-INVENT, we compared
> PROTAC-INVENT
>
> alongside with two widely used docking programs *Glide* \[55\] and
> *Autodock Vina* \[56\] on re-docking of PROTAC into its respective
> PTS. Fourteen published PROTAC PTS crystal structures were used as our
> test set. For *Glide* and *Vina* docking, the initial 3D conformation
> of 14 PROTAC molecules were converted from SMILES by using the
> *LigPrep* module of *Schrodinger* and then docked into the complex
> structure composed by E3 ligase and POI. Initial 3D conformations of
> PROTAC in the workflow of PROTAC-INVENT were optimized by *ROCS*,
> which serves the purpose of bringing the generated PROTAC to the
> vicinity of the reference ligand. As a rule of thumb, initial
> conformation and starting pose of a compound will certainly affect the
> subsequent docking process (such as elapsed time, docking score
> *etc*). For making a fair comparison, we would like to eliminate the
> effect of the initial conformation on docking results. Therefore, in
> addition to using *LigPrep* for conformation generation, we also
> applied*Omega* \[32\] and *ROCS* software \[33\] to generate multiple
> conformations and obtain best aligned 3D conformation as the starting
> conformation of docking.

  -----------------------------------------------------------------------
  *3D based generative    \|                      9
  PROTAC linker design*                           
  ----------------------- ----------------------- -----------------------

  -----------------------------------------------------------------------

> **Table 1.** Comparison on the ratio of compound having desirable
> properties
>
> **Score Component: Average Score**

+-------------+-------------+-------------+-------------+-------------+
|             | **Model**   | **Total     | **Average   | > **Ratio** |
|             |             | unique      | score**     |             |
|             |             | count**     | *\>*        |             |
|             |             |             | **0.8**     |             |
+=============+=============+=============+=============+=============+
| > **Link    | **No        | 17 673      | > 12 441    | > **0.704** |
| > Length**  | C           |             |             |             |
|             | onstraint** |             |             |             |
+-------------+-------------+-------------+-------------+-------------+
|             | **7--9**    | 17 869      | > 9692      | > **0.542** |
+-------------+-------------+-------------+-------------+-------------+
|             | **7--11**   | 17 326      | > 11 574    | > **0.668** |
+-------------+-------------+-------------+-------------+-------------+
|             | **7--13**   | 18 108      | > 12 953    | > **0.715** |
+-------------+-------------+-------------+-------------+-------------+
|             | **7--15**   | 17 714      | > 12 164    | > **0.687** |
+-------------+-------------+-------------+-------------+-------------+
|             | **No RL**   | 5124        | > 2586      | > **0.505** |
+-------------+-------------+-------------+-------------+-------------+

> **Score Component: Docking Score**

+-------------+-------------+-------------+-------------+-------------+
|             | **Model**   | **Total     | > **Docking | > **Ratio** |
|             |             | unique      | > Score**   |             |
|             |             | count**     | > *\<*      |             |
|             |             |             | > −**10.9** |             |
+=============+=============+=============+=============+=============+
| > **Link    | **No        | 17 673      |             | > **0.536** |
| > length**  | c           |             |             |             |
|             | onstraint** |             |             |             |
+-------------+-------------+-------------+-------------+-------------+
|             |             |             | > 9470      |             |
+-------------+-------------+-------------+-------------+-------------+
|             | **7--9**    | 17 869      | > 13 005    | > **0.728** |
+-------------+-------------+-------------+-------------+-------------+
|             | **7--11**   | 17 326      | > 11 831    | > **0.683** |
+-------------+-------------+-------------+-------------+-------------+
|             | **7--13**   | 18 108      | > 11 610    | > **0.641** |
+-------------+-------------+-------------+-------------+-------------+
|             | **7--15**   | 17 714      | > 11 738    | > **0.663** |
+-------------+-------------+-------------+-------------+-------------+
|             | **No RL**   | 5124        | > 1782      | > **0.348** |
+-------------+-------------+-------------+-------------+-------------+

> **Score Component: *PS***

+-------------+-------------+-------------+-------------+-------------+
|             | **Model**   | **Total     | > ***PS**   | **Ratio**   |
|             |             | unique      | > \>*       |             |
|             |             | count**     | > **0.9**   |             |
+=============+=============+=============+=============+=============+
| > **Link    | **No        | 17 673      | > 10 556    | **0.597**   |
| > Length**  | c           |             |             |             |
|             | onstraint** |             |             |             |
+-------------+-------------+-------------+-------------+-------------+
|             | **7--9**    | 17 869      | > 12 891    | **0.721**   |
+-------------+-------------+-------------+-------------+-------------+
|             | **7--11**   | 17 326      | > 12 705    | **0.733**   |
+-------------+-------------+-------------+-------------+-------------+
|             | **7--13**   | 18 108      | > 12 390    | **0.684**   |
+-------------+-------------+-------------+-------------+-------------+
|             | **7--15**   | 17 714      | > 12 175    | **0.687**   |
+-------------+-------------+-------------+-------------+-------------+
|             | **No RL**   | 5124        | > 2554      | **0.498**   |
+-------------+-------------+-------------+-------------+-------------+

> Note:a) '**No RL**' refers to random sampling from the prior model.
> Bold emphasis represents the ratio of the number of molecules meeting
> the scoring component threshold to the total number of generated
> molecules under different constraints.

![](vertopal_6f7e04f02c2f4d5095038d2511b10a73/media/image11.png){width="6.501388888888889in"
height="3.622221128608924in"}

> Figure 10. The RMSDs of re-docking PROTAC for various methods.
>
> The RMSD and elapsed time of re-docking were evaluated and the results
> can be seen in Figure 10 and Table 2. As shown in Figure 10, for
> *Glide* docking using *ROCS* results as starting confor-mation could
> slightly improve performance, in which, for eight
>
> structures, the RMSD of docked conformations of using *ROCS* is better
> than using *LigPrep* alone. While for *Vina*, no advantage of using
> *ROCS* in docking protocol was observed. Overall, best dock-ing
> solution in 12 of 14 structures comes from PROTAC-INVENT,

+-----------------------+-----------------------+-----------------------+
| > 10                  | \|                    | > *Li* et al.         |
+=======================+=======================+=======================+
+-----------------------+-----------------------+-----------------------+

![](vertopal_6f7e04f02c2f4d5095038d2511b10a73/media/image12.png){width="6.5in"
height="7.022222222222222in"}

> Figure 11. The comparison of docking performance of different docking
> protocal in re-docking procedure of X-ray PROTAC of 6HAX. The run time
> and RMSD values to the X-ray conformation were used as assessment
> criteria.
>
> which clearly shows the efficiency of the constrained docking protocol
> used in PROTAC-INVENT. In terms of computation speed (as shown in
> Table 2), the docking protocol in PROTAC-INVENT also out-performed all
> other methods, regardless of using *ROCS* or not in the docking
> protocol. The average elapsed time of our approach is 0.65 minute,
> while that of *Glide* and *Vina* docking in different starting
> conformation could reach 3.46--3.92 minutes and 3.12--3.54 minutes,
> respectively (Table 2).
>
> The 6HAX (crystal structure of BAF \[52\]) system was taken as the
> example (as shown in Figure 11). Among five docking protocols in
> Figure 11, the average running time of *Glide* and *Vina* were 1.85
> and 1.95 minutes, while that of PROTAC-INVENT is only 0.55 minute. The
> RMSD of PROTAC-INVENT is 0.157 Å,
>
> while those of Glide (4.205--4.512 Å) and Vina (0.550--0.666 Å) are
> much higher. These results demonstrate that the constrained docking
> protocol of PROTAC-INVENT is a better choice for doing PROTAC docking
> when it was hypothesized that a PROTAC may mimic a known PROTAC
> binding pose. Our approach performs better in this scenario due to two
> reasons: First, it adopts a linker conformation to make sure two
> terminal warheads mimics those of template as close as possible;
> Second, the 'local only' mode of *Vina* docking was used to fine tune
> the initial conformation so that the docking speed was greatly
> improved while the docking pose of PROTAC molecule can still fit with
> the binding pocket.
>
> In summary, our results demonstrate that the proposed PROTAC-INVENT
> protocol can effectively generate SMILES of

  -----------------------------------------------------------------------
  *3D based generative    \|                      11
  PROTAC linker design*                           
  ----------------------- ----------------------- -----------------------

  -----------------------------------------------------------------------

> **Table 2.** The computation speed for various methods (minute)

+-----------+-----------+-----------+-----------+-----------+-----------+
| > **PDB   | > **      |           | >         |           | >         |
| > code**  | *I*nitial |           | **Initial |           |  **PROTAC |
|           | > ligand  |           | > ligand  |           | -INVENT** |
|           | > poses   |           | > poses   |           |           |
|           | > from    |           | > from    |           |           |
|           | > ROCS**  |           | >         |           |           |
|           |           |           | LigPrep** |           |           |
+===========+===========+===========+===========+===========+===========+
|           | > ***Gli  | ***ADV*** | > ***Gli  | ***ADV*** |           |
|           | de*(SP)** |           | de*(SP)** |           |           |
+-----------+-----------+-----------+-----------+-----------+-----------+
| > **5     | > 2.86    | 3.41      | > 3.07    | 3.36      | > 0.72    |
| > T35**   |           |           |           |           |           |
+-----------+-----------+-----------+-----------+-----------+-----------+
| >         | > 1.07    | 1.42      | > 1.48    | 1.27      | > 0.69    |
|  **6BN7** |           |           |           |           |           |
+-----------+-----------+-----------+-----------+-----------+-----------+
| >         | > 1.16    | 1.43      | > 0.97    | 1.70      | > 0.66    |
|  **6BOY** |           |           |           |           |           |
+-----------+-----------+-----------+-----------+-----------+-----------+
| >         | > 1.92    | 2.34      | > 2.28    | 2.07      | > 0.54    |
|  **6HAY** |           |           |           |           |           |
+-----------+-----------+-----------+-----------+-----------+-----------+
| >         | > 1.86    | 2.15      | > 1.67    | 1.38      | > 0.59    |
|  **6HR2** |           |           |           |           |           |
+-----------+-----------+-----------+-----------+-----------+-----------+
| >         | > 2.61    | 2.44      | > 3.35    | 2.13      | >         |
|  **6W7O** |           |           |           |           |  **0.80** |
+-----------+-----------+-----------+-----------+-----------+-----------+
| >         | >         | **9.89**  | > 8.05    | 9.41      | > 0.68    |
|  **6ZHC** |  **8.34** |           |           |           |           |
+-----------+-----------+-----------+-----------+-----------+-----------+
| >         | > 5.82    | 4.88      | > 5.77    | 4.69      | > 0.70    |
|  **7JTO** |           |           |           |           |           |
+-----------+-----------+-----------+-----------+-----------+-----------+
| >         | > 2.25    | 1.93      | > 2.03    | 1.48      | > 0.61    |
|  **7JTP** |           |           |           |           |           |
+-----------+-----------+-----------+-----------+-----------+-----------+
| >         | > 2.55    | 3.57      | > 2.78    | 2.64      | > 0.60    |
|  **7KHH** |           |           |           |           |           |
+-----------+-----------+-----------+-----------+-----------+-----------+
| >         | > 2.92    | 2.10      | > 3.00    | 1.84      | > 0.61    |
|  **7Q2J** |           |           |           |           |           |
+-----------+-----------+-----------+-----------+-----------+-----------+
| >         | > 5.21    | 6.10      | >         | 3.81      | > 0.65    |
|  **6HMO** |           |           | **10.38** |           |           |
+-----------+-----------+-----------+-----------+-----------+-----------+
| >         | > 7.93    | 5.84      | > 8.27    | 5.99      | > 0.65    |
|  **6W8I** |           |           |           |           |           |
+-----------+-----------+-----------+-----------+-----------+-----------+
| >         | > 1.91    | 2.02      | > 1.78    | 1.89      | > 0.55    |
|  **6HAX** |           |           |           |           |           |
+-----------+-----------+-----------+-----------+-----------+-----------+
| > **      | >         | **3.54**  | >         | **3.12**  | >         |
| Average** |  **3.46** |           |  **3.92** |           |  **0.65** |
+-----------+-----------+-----------+-----------+-----------+-----------+

> Bold emphasis represents the Highest and Average run time under
> different docking protocol.
>
> novel PROTAC molecules along with their hypothesized 3D binding poses
> which mimick that of reference PROTAC. However, as the published
> PROTAC ternary complex structures are still very limited, this could
> be a hurdle for applying this method on

+-----------------------------------------------------------------------+
| > • PROTAC-INVENT can be used as an efficient and robust docking      |
| > protocol for PROTACs.                                               |
+=======================================================================+
+-----------------------------------------------------------------------+

> systems which do not have any PTS. Recently, some efforts \[7,\
> 57--60\] were reported on prediction of PTS. Advancement in this
>
> area could largely help to expand the application domain of our
> method.
>
> CONCLUSIONS
>
> In current work, a novel 3D linker generative model, PROTAC-INVENT,
> was proposed to rationally design PROTAC molecule. Previous analysis
> has shown that linker structure can largely affect the degradation
> efficacy of PROTAC molecules. So far, most of existing linker
> generation method of PROTAC can only design linker on either 1D SMILES
> or 2D graph level, which doesn't take into account the 3D information
> of the E3 ligase/POI complex structure. In PROTAC-INVENT model, the
> protein complex struc-tures were introduced and for the first time,
> the design of PROTAC linker was done at 3D level in which the 3D
> binding conformation PROTAC can be generated and the conformation of
> POI warhead
>
> ABBREVIATIONS\
> PTS, PROTAC ternary structure; ADV , *Autodock Vina*; PS, Post-Dock
> shape similarity of warheads; Warheads, warhead and E3-ligand; RPC,
> re-combined PROTAC conformation; DPC, desired PROTAC conformation. RL,
> Reinforcement Learning.
>
> SUPPLEMENTARY DATA\
> are available online at.
>
> CONFLICT OF INTEREST\
> The authors declare that they have no competing interests.
>
> and E3-ligand can be adjusted to the PTS pocket to mimic those of
> reference structure. The re-docking study of PROTAC-INVENT FUNDING
>
> on known PROTAC crystal structures was also carried out and the
> performance was compared with Glide and *Vina*, our model achieved
> best RMSD and computation speed among these meth-
>
> This work is supported by the fundings from the Basic and Applied
> Basic Research Foundation of Guangzhou (NO. 202102080402), the Pearl
> River Recruitment Program of Talents (NO. 2021CX020227)
>
> ods. and the Key Research and Development Program of Guangdong

Province (NO. 2010A080813002).

+-----------------------------------------------------------------------+
| > **Key Points**                                                      |
| >                                                                     |
| > • PROTAC-INVENT is a novel 3D linker generative model which can not |
| > only generate SMILES of PROTAC but also its 3D putative binding     |
| > conformation coupled with the target protein and the E3 ligase.     |
| >                                                                     |
| > • PROTAC-INVENT is trained jointly with the reinforce-ment learning |
| > approach to bias the generation of PRO-TAC structures toward        |
| > pre-defined 2D and 3D based properties.                             |
+=======================================================================+
+-----------------------------------------------------------------------+

> DATA AVAILABILITY STATEMENT\
> The scripts and model building can be found in the GitHub repository
> ().
>
> REFERENCES\
> 1. Sakamoto KM, Kim KB, Kumagai A, *et al.* Protacs: chimeric
> molecules that target proteins to the Skp1-Cullin-F box complex for
> ubiquitination and degradation. *Proc Natl Acad Sci U S A*
> 2001;**98**:8554--9.

+-----------------+-----------------+-----------------+-----------------+
| > 12            | \|              | > *Li* et al.   |                 |
+=================+=================+=================+=================+
| > 2\. Deshaies  |                 |                 | > 23\. He J,    |
| > RJ. Protein   |                 |                 | > You H,        |
| > degradation:  |                 |                 | > Sandström E,  |
| > prime time    |                 |                 | > *et al.*      |
| > for PROTACs.  |                 |                 | > Molecular     |
| > *Nat Chem     |                 |                 | > optimization  |
| > Biol*         |                 |                 | > by cap-turing |
| > 2015          |                 |                 | > chemist's     |
| ;**11**:634--5. |                 |                 | > intuition     |
|                 |                 |                 | > using deep    |
| 3\. Dale B,     |                 |                 | > neural        |
| Cheng M, Park   |                 |                 | > networks. *J  |
| KS, *et al.*    |                 |                 | > Chem*         |
| Advancing       |                 |                 | > 202           |
| targeted        |                 |                 | 1;**13**:1--17. |
| protein         |                 |                 |                 |
+-----------------+-----------------+-----------------+-----------------+

> degradation for cancer therapy. *Nat Rev Cancer* 2021;**21**:638--54.
> 4. Pettersson M, Crews CM. PROteolysis TArgeting chimeras (PRO-
>
> 24\. Arús-Pous J, Blaschke T, Ulander S, *et al.* Exploring the GDB-13
> chemical space using deep generative models. *J Chem* 2019;**11**:
>
> TACs) --- past, present and future. *Drug Discov Today Technol* 1--14.
>
> 2019;**31**:15--27. 25. Yang Y, Zheng S, Su S, *et al.* SyntaLinker:
> automatic fragment
>
> 5\. Lai AC, Crews CM. Induced protein degradation: an emerging drug
> discovery paradigm. *Nat Rev Drug Discov* 2017;**16**:101--14. 6. Bai
> L, Zhou H, Xu R, *et al.* A potent and selective small-molecule
> degrader of STAT3 achieves complete tumor regression in vivo.
>
> linking with deep conditional transformer neural networks. *Chem Sci*
> 2020;**11**:8312--22.
>
> 26\. Zheng S, Tan Y, Wang Z, *et al.* Accelerated rational PROTAC
> design via deep learning and molecular simulations. *Nat Mach*
>
> *Cancer Cell* 2019;**36**:498--511.e17. *Intell* 2022;**4**:739--48.
>
> 7\. Tunjic TM, Weber N, Brunsteiner M. Computer aided drug design in
> the development of proteolysis targeting chimeras. *Comput Struct
> Biotechnol J* 2023;**21**:2058--67.
>
> 8\. Purohit R, Rajasekaran R, Sudandiradoss C, *et al.* Studies on
> flex-
>
> 27\. Imrie F, Bradley AR, Van Der Schaar M, *et al.* Deep generative
> models for 3D linker design. *J Chem Inf Model* 2020;**60**:1983--95.
> 28. Igashov I, Stärk H, Vignac C, *et al.* Equivariant 3D-conditional
> dif- fusion models for molecular linker design. 2022. arXiv preprint
>
> ibility and binding affinity of Asp25 of HIV-1 protease mutants.
> arXiv:2210.05274.
>
> *Int J Biol Macromol* 2008;**42**:386--91.
>
> 9\. Kumar A, Rajendran V, Sethumadhavan R, *et al.* In silico pre-
>
> 29\. Guo J, Knuth F, Margreitter C, *et al.* Link-INVENT: genera-tive
> linker design with reinforcement learning. *Digital Discovery*
>
> diction of a disease-associated STIL mutant and its affect on
> 2023;**2**:392--408.
>
> the recruitment of centromere protein J (CENPJ). *FEBS Open Bio* 30.
> Lipinski CA, Lombardo F, Dominy BW, *et al.* Experimental
>
> 2012;**2**:285--93. and computational approaches to estimate
> solubility and
>
> 10\. Tanwar G, Mazumder AG, Bhardwaj V, *et al.* Target
> iden-tification, screening and in vivo evaluation of pyrrolone-fused
> benzosuberene compounds against human epilepsy using zebrafish model
> of pentylenetetrazol-induced seizures. *Sci Rep*
>
> permeability in drug discovery and development settings. *Adv Drug
> Deliv Rev* 2001;**46**:3--26.
>
> 31\. Blaschke T, Arús-Pous J, Chen H, *et al.* REINVENT 2.0: an AI
> tool for de novo drug design. *J Chem Inf Model* 2020;**60**:5918--22.
>
> 2019;**9**:7904. 32. Hawkins PCD, Skillman AG, Warren GL, *et al.*
> Conformer gener-
>
> 11\. Bhardwaj V, Purohit R. Computational investigation on effect of
> mutations in PCNA resulting in structural perturbations and inhibition
> of mismatch repair pathway. *J Biomol Struct Dyn*
> 2020;**38**:1963--74.
>
> 12\. Kumar Bhardwaj V, Purohit R, Kumar S. Himalayan bioactive
>
> ation with OMEGA: algorithm and validation using high quality
> structures from the protein databank and Cambridge structural
> database. *J Chem Inf Model* 2010;**50**:572--84.
>
> 33\. Hawkins PCD, Skillman AG, Nicholls A. Comparison of
> shape-matching and docking as virtual screening tools. *J Med Chem*
>
> molecules as potential entry inhibitors for the human immun-
> 2007;**50**:74--82.
>
> odeficiency virus. *Food Chem* 2021;**347**:128932. 34. Schrödinger
> Release 2022-4. *MacroModel, Schrödinger*. New York,
>
> 13\. Bhardwaj VK, Oakley A, Purohit R. Mechanistic behavior and NY:
> LLC, 2021.
>
> subtle key events during DNA clamp opening and closing in T4
> bacteriophage. *Int J Biol Macromol* 2022;**208**:11--9.
>
> 14\. Liu Z, Hu X, Wang Q, *et al.* Design and synthesis of EZH2-based
> PROTACs to degrade the PRC2 complex for target-ing the noncatalytic
> activity of EZH2. *J Med Chem* 2021;**64**:
>
> 35\. Guo J, Janet JP, Bauer MR, *et al.* DockStream: a docking wrapper
> to enhance de novo molecular design. *J Chem* 2021;**13**:1--21.
>
> 36\. Tang S, Chen R, Lin M, *et al.* Accelerating AutoDock VINA with
> GPUs. *ChemRxiv* 2022;**27**:3041.
>
> 37\. Huey R, Morris GM, Forli S. Using AutoDock 4 and AutoDock

+--------+--------+--------+--------+--------+--------+--------+--------+
| > 282  | Vina   | with   | Au     | a      | tut    | 2012;  | h      |
| 9--48. |        |        | toDock |        | orial. | 1--32. | ttp:// |
|        |        |        | Tools: |        |        |        |        |
+========+========+========+========+========+========+========+========+
+--------+--------+--------+--------+--------+--------+--------+--------+

> 15\. Weng G, Cai X, Cao D, *et al.* PROTAC-DB 2.0: an updated database
> of PROTACs. *Nucleic Acids Res* 2023;**51**:D1367--72.
>
> 16\. Smith BE, Wang SL, Jaime-Figueroa S, *et al.* Differential PROTAC
> substrate specificity dictated by orientation of recruited E3 lig-ase.
> *Nat Commun* 2019;**10**:1--13.
>
> .
>
> 38\. Morris GM, Ruth H, Lindstrom W, *et al.* Software news and
> updates AutoDock4 and AutoDockTools4: automated dock-ing with
> selective receptor flexibility. *J Comput Chem* 2009;**30**:
>
> 17\. Bemis TA, La CJJ, Burkart MD. Unraveling the role of linker
> design 2785--91.
>
> in proteolysis targeting chimeras. *J Med Chem* 2021;**64**:8042--52.
> 18. Mohler ML, Sikdar A, Ponnusamy S, *et al.* An overview of next-
> generation androgen receptor-targeted therapeutics in devel- opment
> for the treatment of prostate cancer. *Int J Mol Sci*
>
> 39\. Pike A, Williamson B, Harlfinger S, *et al.* Optimising
> proteolysis-targeting chimeras (PROTACs) for oral drug delivery: a
> drug metabolism and pharmacokinetics perspective. *Drug Discov Today*
> 2020;**25**:1793--800.
>
> 2021;**22**:2124. 40. Cantrill C, Chaturvedi P, Rynn C, *et al.*
> Fundamental aspects of
>
> 19\. Li Y, Hu J, Wang Y, *et al.* DeepScaffold: a comprehensive tool
> DMPK optimization of targeted protein degraders. *Drug Discov*
>
> for scaffold-based de novo drug discovery using deep learning. *J*
> *Today* 2020;**25**:969--82.
>
> *Chem Inf Model* 2020;**60**:77--91.
>
> 20\. Li Y, Zhang L, Liu Z. Multi-objective de novo drug design with
>
> 41\. Schiemer J, Horst R, Meng Y, *et al.* Snapshots and ensembles of
> BTK and cIAP1 protein degrader ternary complexes. *Nat Chem*
>
> conditional graph generative model. *J Chem* 2018;**10**:33. *Biol*
> 2021;**17**:152--60.
>
> 21\. Skalic M, Jiménez J, Sabbadin D, *et al.* Shape-based generative
> 42. Cummins DJ, Bell MA. Integrating everything: the molecule modeling
> for de novo drug design. *J Chem Inf Model* 2019;**59**: selection
> toolkit, a system for compound prioritization in drug 1205--14.
> discovery. *J Med Chem* 2016;**59**:6999--7010.
>
> 22\. Sturm N, Sun J, Vandriessche Y, *et al.* Application of
> bioactivity 43. Blaschke T, Engkvist O, Bajorath J, *et al.*
> Memory-assisted rein- profile-based fingerprints for building machine
> learning models. forcement learning for diverse molecular de novo
> design. *J Chem* *J Chem Inf Model* 2019;**59**:962--72.
> 2020;**12**:1--17.

+-----------------+-----------------+-----------------+-----------------+
|                 | *3D based       | \|              | 13              |
|                 | generative      |                 |                 |
|                 | PROTAC linker   |                 |                 |
|                 | design*         |                 |                 |
+=================+=================+=================+=================+
| > 44\. Hendriks | > 52\. Farnaby  |                 |                 |
| > RW, Yuvaraj   | > W, Koegl M,   |                 |                 |
| > S, Kil LP.    | > Roy MJ, *et   |                 |                 |
| > Targeting     | > al.* BAF      |                 |                 |
| > Bruton's      | > complex       |                 |                 |
| > tyro-sine     | >               |                 |                 |
| > kinase in B   | vulnerabilities |                 |                 |
| > cell          | > in cancer     |                 |                 |
| > malignancies. | > demonstrated  |                 |                 |
| > *Nat Rev      | > via           |                 |                 |
| > Cancer*       | >               |                 |                 |
| > 2014;**14**:  | structure-based |                 |                 |
| > 219--32.      | > PROTAC        |                 |                 |
| >               | > design. *Nat  |                 |                 |
| > 45\. Dunleavy | > Chem Biol*    |                 |                 |
| > K, Erdmann T, | > 2019;         |                 |                 |
| > Lenz G.       | **15**:672--80. |                 |                 |
| > Targeting the | >               |                 |                 |
| > B-cell        | > 53\. Gadd MS, |                 |                 |
| > receptor      | > Testa A,      |                 |                 |
| > pathway in    | > Lucas X, *et  |                 |                 |
| > diffuse large | > al.*          |                 |                 |
| > B-cell        | > Structural    |                 |                 |
| > lymphoma.     | > basis of      |                 |                 |
| > *Cancer Treat | > PROTAC        |                 |                 |
| > Rev*          | > cooperative   |                 |                 |
| > 201           | > recognition   |                 |                 |
| 8;**65**:41--6. | > for selective |                 |                 |
| >               | > protein       |                 |                 |
| > 46\. Byrd JC, | > degradation.  |                 |                 |
| > Furman RR,    | > *Nat Chem     |                 |                 |
| > Coutre SE,    | > Biol*         |                 |                 |
| > *et al.*      | > 2017;         |                 |                 |
| > Targeting BTK | **13**:514--21. |                 |                 |
| > with          | >               |                 |                 |
| > Ibrutinib in  | > 54\. Bowers   |                 |                 |
| > relapsed      | > KJ, Sacerdoti |                 |                 |
| > chronic       | > FD, Salmon    |                 |                 |
| > lymphocytic   | > JK, *et al.*  |                 |                 |
| > leukemia. *N  | > Molecular     |                 |                 |
| > Engl J Med*   | > dyna          |                 |                 |
| > 2013;**369**: | mics---Scalable |                 |                 |
| > 32--42.       | > algorithms    |                 |                 |
|                 | > for molecular |                 |                 |
|                 | > dynamics      |                 |                 |
|                 | > simulations   |                 |                 |
|                 | > on commodity  |                 |                 |
|                 | > clusters.     |                 |                 |
|                 | > *Proceedings  |                 |                 |
|                 | > of the 2006   |                 |                 |
|                 | > ACM/IEEE      |                 |                 |
|                 | > conference*   |                 |                 |
+-----------------+-----------------+-----------------+-----------------+

> 47\. Wang ML, Rule S, Martin P, *et al.* Targeting BTK with Ibrutinib
> in relapsed or refractory mantle-cell lymphoma. *N Engl J Med*
> 2013;**369**:507--16.
>
> *on Supercomputing -- SC'06* 2006; 84.
>
> 55\. Friesner RA, Banks JL, Murphy RB, *et al.* Glide: a new approach
> for rapid, accurate docking and scoring. 1. Method and assessment
>
> 48\. Qiu H, Liu-Bujalski L, Caldwell RD, *et* *al.* Discovery of of
> docking accuracy. *J Med Chem* 2004;**47**:1739--49.
>
> potent, highly selective covalent irreversible BTK inhibitors 56.
> Trott O, Olson AJ. AutoDock Vina: improving the speed and
>
> from a fragment hit. *Bioorg Med* *Chem* *Lett* 2018;**28**: accuracy
> of docking with a new scoring function, efficient opti-
>
> 2939--44. mization, and multithreading. *J Comput Chem*
> 2009;**29**:NA-NA.
>
> 49\. Pan Z, Scheerens H, Li S-J, *et al.* Discovery of selective
> irreversible 57. Drummond ML, Williams CI. In silico modeling of
> PROTAC-
>
> inhibitors for Bruton's tyrosine kinase. *Chem Med Chem* 2007;**2**:
> mediated ternary complexes: validation and application. *J Chem*
>
> 58--61. *Inf Model* 2019;**59**:1634--44.
>
> 50\. Wu H, Huang Q, Qi Z, *et al.* Irreversible inhibition of BTK
> kinase 58. Rao A, Tunjic TM, Brunsteiner M, *et al.* Bayesian
> optimization
>
> by a novel highly selective inhibitor CHMFL-BTK-11 suppresses for
> ternary complex prediction (BOTCP). *Artific Intell Life Sci*
>
> inflammatory response in rheumatoid arthritis model. *Sci Rep*
> 2023;**3**:100072.
>
> 2017;**7**:1--10. 59. Gao H, Sun X, Rao Y. PROTAC technology:
> opportunities and
>
> 51\. Huang H-T, Dobrovolsky D, Paulk J, *et al.* A chemoproteomic
> approach to query the degradable kinome using a multi-kinase degrader.
> *Cell Chem Biol* 2018;**25**:88--99.e6.
>
> challenges. *ACS Med Chem Lett* 2020;**11**:237--40.
>
> 60\. Weng G, Li D, Kang Y, *et al.* Integrative modeling of
> PROTAC-mediated ternary complexes. *J Med Chem* 2021;**64**:16271--81.
