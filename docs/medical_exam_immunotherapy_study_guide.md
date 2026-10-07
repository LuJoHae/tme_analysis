# Comprehensive Medical Exam Study Guide: Cancer Immunotherapy & Medication Master Catalog

> **Target Examinations**: USMLE Step 1 / Step 2 CK / Step 3, ABIM Oncology Subspecialty Board Exam, Medical Oncology Specialty Certificate Exam (SCE / MRCP), AMC, PLAB, and Medical Pharmacology Finals.  
> **Key Focus**: Drug mechanisms, receptor-ligand interactions, generic vs. brand names, tumor-specific and agnostic indications, toxicities (irAEs, CRS, ICANS), and clinical management algorithms.

---

## Quick-Reference Nomenclature & Suffix Decoder

Monoclonal antibodies and biological therapies follow standardized International Nonproprietary Name (INN) stem and suffix conventions:

| Substem / Suffix | Meaning & Source Species | Human Sequence % | Immunogenicity Risk | High-Yield Clinical Examples |
| :--- | :--- | :---: | :---: | :--- |
| **`-omab`** | Murine (100% mouse protein) | 0% | Highest | Ibritumomab tiuxetan, historic muromonab-CD3 |
| **`-ximab`** | Chimeric (mouse variable region + human constant region) | ~65% | Moderate | **Rituximab**, Dinutuximab, Infliximab, Cetuximab |
| **`-zumab`** | Humanized (mouse complementarity-determining regions [CDRs] grafted onto human framework) | >90% | Low | **Pembrolizumab**, **Atezolizumab**, Dostarlimab, Trastuzumab |
| **`-umab`** | Fully human (derived from transgenic mice or phage display) | 100% | Lowest | **Nivolumab**, **Ipilimumab**, **Durvalumab**, Cemiplimab |
| **`-cel`** | Cellular immunotherapy (autologous or allogeneic living cells) | N/A | Variable | Tisagenlecleu**cel**, Axicabtagene ciloleu**cel**, Lifileu**cel** |
| **`-vec`** | Viral gene-delivery vector / oncolytic virus | N/A | High | Talimogene laherparep**vec** (T-VEC), Nadofaragene firadeno**vec** |
| **`-cept`** | Receptor-Fc fusion protein | Variable | Low | Eftilagimod alfa (soluble LAG-3), Aflibercept |

---

## 1. Immune Checkpoint Inhibitors (ICIs)

Under physiological conditions, immune checkpoints maintain self-tolerance and limit collateral tissue damage during inflammation. Malignant cells co-opt these inhibitory pathways to evade immune surveillance. Checkpoint blockade releases the "molecular brakes" on tumor-reactive T cells.

```
       ANTIGEN PRESENTATION & PRIMING (Lymph Node)                  PERIPHERAL EFFECTOR PHASE (Tumor Bed)
       
       Dendritic Cell              T Cell                           Tumor Cell                 Effector CD8+ T Cell
      +---------------+      +-----------------+                  +---------------+      +----------------------+
      |   MHC I / II  | ===> |     TCR         |                  |  MHC I - Ag   | ===> |        TCR           |
      |   B7-1 / B7-2 | -X-> | CTLA-4 (OFF)    |                  |    PD-L1      | -X-> |      PD-1 (OFF)      |
      | (CD80 / CD86) |      |   [Ipilimumab]  |                  | (CD274 / B7-H1)|     | [Pembrolizumab/Nivo] |
      +---------------+      +-----------------+                  +---------------+      +----------------------+
                             |  CD28 (ON)      |                                         |  4-1BB / OX40 (ON)   |
                             +-----------------+                                         +----------------------+
```

### A. Anti-PD-1 Monoclonal Antibodies (Receptor Blockade)
* **Target & Mechanism**: Binds **PD-1 (CD279)** on activated T cells, B cells, and NK cells, preventing engagement with its immunosuppressive ligands **PD-L1 (CD274, B7-H1)** and **PD-L2 (CD273, B7-DC)**. Reverses T-cell exhaustion and restores cytotoxic CD8+ T cell tumor killing.

| Generic Name (INN) | Trade / Brand Name | Antibody Type | High-Yield Approved Clinical Indications | Medical Exam Pearls & Buzzwords |
| :--- | :--- | :--- | :--- | :--- |
| **Pembrolizumab** | **Keytruda** | Humanized IgG4 kappa | • **Metastatic Melanoma**<br>• **Non-Small Cell Lung Cancer (NSCLC)** (1st-line monotherapy if PD-L1 TPS $\ge 50\%$; or + chemo regardless of PD-L1)<br>• **MSI-H / dMMR Solid Tumors** (FDA tumor-agnostic approval)<br>• **TMB-High Solid Tumors** ($\ge 10\text{ mut/Mb}$)<br>• Triple-Negative Breast Cancer (TNBC)<br>• Urothelial / Bladder Carcinoma<br>• Renal Cell Carcinoma (RCC, + axitinib or lenvatinib)<br>• Classical Hodgkin Lymphoma (cHL)<br>• Head & Neck Squamous Cell Carcinoma (HNSCC)<br>• Cervical, Esophageal, & Gastric Cancers | • **Tumor-Agnostic Paradigm**: First drug approved based on a genomic biomarker (MSI-H/dMMR) rather than anatomical site.<br>• Requires PD-L1 testing in many 1st-line indications (TPS = Tumor Proportion Score; CPS = Combined Positive Score). |
| **Nivolumab** | **Opdivo** | Fully human IgG4 | • **Melanoma** (adjuvant and metastatic; monotherapy or combined with ipilimumab)<br>• **NSCLC** (alone or + ipilimumab $\pm$ chemo)<br>• **Renal Cell Carcinoma** (+ ipilimumab or cabozantinib)<br>• **Classical Hodgkin Lymphoma** (relapsed post-autologous stem cell transplant & brentuximab)<br>• **HNSCC** (platinum-refractory)<br>• **Urothelial Bladder Cancer**<br>• **Colorectal Cancer (MSI-H / dMMR)** (+ ipilimumab)<br>• **Hepatocellular Carcinoma (HCC)** (+ ipilimumab)<br>• **Esophageal Squamous Cell Carcinoma (ESCC)** | • Combined with Ipilimumab in the landmark **CheckMate-067** trial (dramatically increased survival in metastatic melanoma, but $>50\%$ Grade 3/4 irAEs).<br>• IgG4 backbone minimizes unwanted antibody-dependent cellular cytotoxicity (ADCC) against host T cells. |
| **Cemiplimab** | **Libtayo** | Fully human IgG4 | • **Cutaneous Squamous Cell Carcinoma (CSCC)** (locally advanced or metastatic, not candidates for curative surgery/radiation)<br>• **Basal Cell Carcinoma (BCC)** (previously treated with Hedgehog inhibitor)<br>• **NSCLC** (monotherapy with PD-L1 $\ge 50\%$) | • **Exam Board Favorite**: Standard of care for non-resectable advanced CSCC, which historically had minimal effective systemic therapy options. |
| **Dostarlimab** | **Jemperli** | Humanized IgG4 | • **Mismatch Repair Deficient (dMMR) Recurrent/Advanced Endometrial Cancer**<br>• **dMMR Locally Advanced Rectal Cancer** (neoadjuvant monotherapy) | • **Landmark Trial Result (Cercek et al., NEJM 2022)**: Achieved a **100% complete clinical response rate** without surgery, radiation, or chemotherapy in stage II/III dMMR rectal cancer. |
| **Tislelizumab** | **Tevimbra** | Humanized IgG4 | • **Esophageal Squamous Cell Carcinoma (ESCC)** (unresectable/metastatic after prior systemic chemo)<br>• Gastric and Gastroesophageal Junction (GEJ) adenocarcinoma | • Engineered with an altered Fc region to **abolish Fc$\gamma$R binding on macrophages**, eliminating macrophage-mediated clearance of anti-tumor T cells. |
| **Toripalimab** | **Loqtorzi** | Humanized IgG4 | • **Nasopharyngeal Carcinoma (NPC)** (1st-line with gemcitabine/cisplatin, and recurrent/metastatic post-chemo) | • First FDA-approved drug specifically indicated for nasopharyngeal carcinoma. |
| **Retifanlimab** | **Zynyz** | Humanized IgG4 | • **Merkel Cell Carcinoma (MCC)** (metastatic or recurrent locally advanced) | • Competes with Avelumab and Pembrolizumab in neuroendocrine skin cancers. |
| *Sintilimab (Tyvyt)*<br>*Camrelizumab (Airuika)*<br>*Serplulimab (Hansizhuang)* | *(Various)* | Humanized IgG4 | • Predominantly approved in East Asia (NMPA) for NSCLC, ESCC, HCC, and lymphoma; advancing in global phase III trials. | • High clinical relevance in global clinical trials and international oncology board questions. |

---

### B. Anti-PD-L1 Monoclonal Antibodies (Ligand Blockade)
* **Target & Mechanism**: Binds **PD-L1 (CD274, B7-H1)** expressed on tumor cells and tumor-infiltrating myeloid cells. Blocks interaction with PD-1 and CD80 (B7-1), while leaving the **PD-1 / PD-L2** interaction intact (theoretically leading to fewer autoimmune pneumonitis events than anti-PD-1, though clinically comparable).

| Generic Name (INN) | Trade / Brand Name | Antibody Type | High-Yield Approved Clinical Indications | Medical Exam Pearls & Buzzwords |
| :--- | :--- | :--- | :--- | :--- |
| **Atezolizumab** | **Tecentriq** | Engineered Humanized IgG1 | • **Extensive-Stage Small Cell Lung Cancer (ES-SCLC)** (1st-line with carboplatin + etoposide [IMpower133 trial])<br>• **Hepatocellular Carcinoma (HCC)** (1st-line with Bevacizumab [IMbrave150 trial])<br>• **NSCLC** (adjuvant post-resection and 1st-line metastatic)<br>• Metastatic Urothelial Carcinoma | • Engineered with an **aglycosylated Fc domain (N298A mutation)** to completely abrogate Fc effector function and prevent antibody-dependent destruction of PD-L1+ activated immune cells.<br>• **Atezolizumab + Bevacizumab** replaced sorafenib as the 1st-line standard of care in advanced unresectable HCC. |
| **Durvalumab** | **Imfinzi** | Fully human IgG1 kappa | • **Unresectable Stage III NSCLC** consolidation after definitive concurrent chemoradiotherapy (**PACIFIC Trial**)<br>• **Extensive-Stage SCLC** (with etoposide + platinum)<br>• **Biliary Tract Cancer (BTC)** (cholangiocarcinoma + gallbladder cancer, with gemcitabine/cisplatin [TOPAZ-1 trial])<br>• **Hepatocellular Carcinoma** (combined with Tremelimumab in the STRIDE regimen) | • **The "PACIFIC Regimen"**: Classic board question—12 months of durvalumab consolidation following concurrent chemoradiotherapy in inoperable stage III NSCLC significantly prolongs overall survival. |
| **Avelumab** | **Bavencio** | Fully human IgG1 | • **Advanced Urothelial Carcinoma** (1st-line switch maintenance therapy post-platinum [JAVELIN Bladder 100 trial])<br>• **Merkel Cell Carcinoma (MCC)** (metastatic)<br>• **Advanced Renal Cell Carcinoma** (combined with Axitinib) | • **Unique Mechanism**: Unlike atezolizumab and durvalumab, avelumab retains a **fully functional, competent IgG1 Fc domain**, capable of inducing **NK-cell mediated ADCC** directly against PD-L1+ tumor cells. |
| **Envafolimab** | **Enweida** | Single-domain antibody (Nanobody) Fc fusion | • **MSI-H / dMMR Advanced Solid Tumors** (colorectal, gastric, endometrial) | • First **subcutaneously administered** PD-L1 inhibitor worldwide (camelid heavy-chain single-domain antibody linked to human IgG1 Fc). Avoids IV infusion time. |
| **Sugemalimab** | **Cejemly** | Fully human IgG4 | • Stage III & IV NSCLC, Relapsed/Refractory Extranodal NK/T-cell lymphoma | • Extensively studied in Asian and European cohorts. |
| **Adebrelimab** | **Airuike** | Humanized IgG4 | • Extensive-stage small cell lung cancer (ES-SCLC) | • CAPSTONE-1 phase III trial in 1st-line ES-SCLC. |

---

### C. Anti-CTLA-4 Monoclonal Antibodies
* **Target & Mechanism**: Targets **Cytotoxic T-Lymphocyte-Associated Protein 4 (CTLA-4 / CD152)**. CTLA-4 competes with the costimulatory receptor **CD28** for binding to **B7-1 (CD80)** and **B7-2 (CD86)** on antigen-presenting cells (APCs). CTLA-4 has orders of magnitude higher affinity for B7 than CD28, delivering an inhibitory signal that arrests T-cell activation during early priming in secondary lymphoid organs. CTLA-4 blockade promotes T-cell priming, clonal expansion, and selectively depletes intratumoral **FOXP3+ regulatory T cells (Tregs)** via Fc-mediated ADCC/ADCP.

| Generic Name (INN) | Trade / Brand Name | Target / Class | Approved Clinical Indications | High-Yield Board Exam Points |
| :--- | :--- | :--- | :--- | :--- |
| **Ipilimumab** | **Yervoy** | Fully human IgG1 against CTLA-4 | • **Metastatic Melanoma** (adjuvant and metastatic; monotherapy or + nivolumab)<br>• **Renal Cell Carcinoma (RCC)** (intermediate/poor risk, + nivolumab)<br>• **NSCLC** (1st-line metastatic, + nivolumab $\pm$ 2 cycles histology-based chemo)<br>• **Malignant Pleural Mesothelioma** (1st-line unresectable, + nivolumab)<br>• **Colorectal Cancer (MSI-H / dMMR)** (+ nivolumab)<br>• **Hepatocellular Carcinoma (HCC)** (+ nivolumab)<br>• **Esophageal Squamous Cell Carcinoma** (+ nivolumab) | • **First checkpoint inhibitor approved** by FDA (2011, revolutionized melanoma).<br>• Acts predominantly at the **lymph node priming phase**, whereas anti-PD-1 acts at the **tumor microenvironment effector phase**.<br>• **Classic irAE**: **Autoimmune Hypophysitis (Pituitaritis)**—headaches, visual changes, panhypopituitarism (low ACTH, TSH, LH/FSH) with pituitary enlargement on MRI. |
| **Tremelimumab** | **Imjudo** | Fully human IgG2 against CTLA-4 | • **Unresectable Hepatocellular Carcinoma** (in the **STRIDE regimen**: Single Tremelimumab Regular Interval Durvalumab)<br>• **Metastatic NSCLC** (combined with durvalumab + platinum chemotherapy) | • **STRIDE concept**: A single priming dose of tremelimumab ($300\text{ mg}$) is given once upfront to stimulate T-cell repertoire diversity, followed by continuous durvalumab monotherapy, maximizing efficacy while minimizing cumulative CTLA-4 autoimmune toxicity. |

---

### D. Anti-LAG-3, TIGIT, TIM-3 & Novel/Bispecific Checkpoint Blockade

| Generic Name | Target Receptor & Mechanism | Clinical Status & Brand Name | Indications / Key Trials | High-Yield Exam Context |
| :--- | :--- | :--- | :--- | :--- |
| **Relatlimab** | **LAG-3 (CD223)**: Binds MHC Class II on APCs with higher affinity than CD4, suppressing T-cell activation. Often co-expressed with PD-1 on severely exhausted CD8+ T cells. | **Opdualag** (Fixed-dose combination: **Relatlimab 160 mg + Nivolumab 480 mg**) | • **Unresectable or Metastatic Melanoma** (1st-line, RELATIVITY-047 trial) | • First-in-class LAG-3 inhibitor.<br>• Fixed-dose co-infusion provides superior progression-free survival (PFS) compared to nivolumab monotherapy with significantly less toxicity than nivolumab + ipilimumab. |
| **Tiragolumab** | **TIGIT**: Competes with CD226 for binding to CD155 (PVR) and CD112, suppressing CD8+ T and NK cell cytotoxicity while promoting Treg activity. | Investigational (Phase III: SKYSCRAPER trials) | • NSCLC, ES-SCLC, ESCC (combined with atezolizumab) | • TIGIT inhibition is synergistic with PD-L1 blockade; trials highlight the dual role of restoring effector function and destabilizing Tregs. |
| **Cadonilimab** | **Dual Checkpoint: Bispecific Anti-PD-1 / Anti-CTLA-4** | **Kaitani** (Approved in China / Phase III global) | • Relapsed/metastatic **Cervical Cancer**, Gastric/GEJ adenocarcinoma | • First approved bispecific ICI worldwide. Tetravalent design exhibits higher avidity in tumor tissue with high density of both PD-1 and CTLA-4, reducing off-target systemic irAEs. |
| **Ivonescimab** | **Bispecific Anti-PD-1 / Anti-VEGF-A** | Approved (HARMONi trial) | • **EGFR-mutated NSCLC** progressing after EGFR-TKI; 1st-line NSCLC | • Combines immune checkpoint disinhibition with anti-angiogenic vascular normalization, reversing the immunosuppressive tumor endothelium. |

---

## 2. Bispecific T-Cell Engagers (BiTEs) & Redirecting Antibodies

Bispecific antibodies possess two distinct antigen-binding Fab arms: one arm binds a tumor-associated surface antigen, while the other binds **CD3$\epsilon$** on the T-cell receptor (TCR) complex. This physically forces a cytotoxic synapse between any polyclonal T cell and the tumor cell **independently of MHC Class I presentation or TCR specificity**, triggering perforin/granzyme-mediated apoptosis.

```
                         [ BISPECIFIC ENGAGER ]
                       Fab (Tumor)     Fab (T Cell)
                           \             /
                            \           /
                             [-- Fc --]
                             /         \
                            v           v
                     +------------+   +------------+
                     | Tumor Cell |   | CD8+ T Cell|
                     |  Antigen   |   |    CD3     |
                     +------------+   +------------+
                            ===> CYTOTOXIC SYNAPSE <===
                       Perforin / Granzyme B Release -> Lysis
```

| Generic Name (INN) | Trade / Brand Name | Molecular Targets | Approved Clinical Indications | Classic Exam Hallmarks & Side Effects |
| :--- | :--- | :--- | :--- | :--- |
| **Blinatumomab** | **Blincyto** | **CD19 $\times$ CD3** (Pure BiTE; lacks Fc domain; short $t_{1/2} \approx 2\text{ hours}$) | • **B-cell Precursor Acute Lymphoblastic Leukemia (B-ALL)** in hematologic remission with **Minimal Residual Disease (MRD) $\ge 0.1\%$**<br>• Relapsed or refractory B-ALL | • Administered as a **continuous IV infusion over 28 days** via portable pump due to its ultra-short half-life.<br>• **Black Box Warnings**: Cytokine Release Syndrome (CRS) and severe neurological toxicities (ICANS). |
| **Tebentafusp** | **Kimmtrak** | **gp100 peptide-HLA-A\*02:01 $\times$ CD3** (ImmTAC: soluble TCR fused to anti-CD3 scFv) | • **Unresectable or Metastatic Uveal Melanoma** in **HLA-A\*02:01-positive adults** | • First TCR-based bispecific therapeutic approved.<br>• **Board Exam Testing Requirement**: Strict requirement for **HLA-A\*02:01 typing** before prescribing; will not work if the patient does not express this specific MHC allele.<br>• Classic toxicity: Cytokine-mediated erythema, edema, and hypotension. |
| **Tarlatamab** | **Imdelltra** | **DLL3 (Delta-like ligand 3) $\times$ CD3** (Half-life extended with Fc domain) | • **Extensive-Stage Small Cell Lung Cancer (ES-SCLC)** with progression on or after platinum chemo | • DLL3 is selectively expressed on neuroendocrine tumor cells (SCLC) with minimal expression in normal tissues.<br>• Major landmark: First BiTE approved for a major common solid tumor indication. |
| **Teclistamab** | **Tecvayli** | **BCMA (B-Cell Maturation Antigen) $\times$ CD3** | • **Relapsed or Refractory Multiple Myeloma (RRMM)** after $\ge 4$ prior lines (including IMiD, PI, and anti-CD38) | • High risk of severe opportunistic infections (Pneumocystis, CMV, fungal) due to complete eradication of normal plasma cells; requires intravenous immunoglobulin (IVIG) and antimicrobial prophylaxis. |
| **Elranatamab** | **Elrexfio** | **BCMA $\times$ CD3** | • **Relapsed or Refractory Multiple Myeloma** after $\ge 4$ prior lines | • Administered subcutaneously with a step-up dosing schedule to mitigate the severity of Cytokine Release Syndrome. |
| **Talquetamab** | **Talvey** | **GPRC5D (G-protein coupled receptor family C group 5 member D) $\times$ CD3** | • **Relapsed or Refractory Multiple Myeloma** | • **High-Yield Unique Toxicities**: GPRC5D is heavily expressed in keratinized tissues (skin, nail bed, filiform papillae of tongue). Patients develop **severe dysgeusia (loss of taste / metallic taste)** leading to anorexia/weight loss, dry mouth, skin rash, and **nail dystrophy / onycholysis**. |
| **Mosunetuzumab** | **Lunsumio** | **CD20 $\times$ CD3** (Full-length humanized IgG1 with silenced Fc) | • **Relapsed or Refractory Follicular Lymphoma** after $\ge 2$ lines of systemic therapy | • "Off-the-shelf" T-cell redirecting therapy available immediately without the vein-to-vein manufacturing delay associated with CAR-T cells. |
| **Glofitamab** | **Columvi** | **CD20 $\times$ CD3** (Novel **2:1 bivalent CD20** to 1 CD3 format) | • **Diffuse Large B-Cell Lymphoma (DLBCL)** or Large B-Cell Lymphoma after $\ge 2$ lines | • **2:1 structure**: Contains two Fab regions binding CD20 and one binding CD3, providing ultra-high avidity for CD20.<br>• Requires mandatory pre-treatment with **Obinutuzumab** (anti-CD20) to deplete circulating B cells and mitigate CRS. |
| **Epcoritamab** | **Epkinly** | **CD20 $\times$ CD3** | • **DLBCL** and **Follicular Lymphoma** | • Delivered as a **subcutaneous injection**, reducing infusion-related hypersensitivity and bed occupancy. |

---

## 3. Cellular Immunotherapies (Adoptive Cell Transfer)

Cellular therapies involve harvesting, genetically engineering (or expanding), and infusing immune effector cells into the patient. Most protocols mandate **lymphodepleting conditioning chemotherapy** (typically Fludarabine + Cyclophosphamide) prior to infusion to eliminate immunosuppressive Tregs and host lymphocytes, clearing "cytokine sinks" (elevating homeostatic IL-7 and IL-15).

```
         AUTOLOGOUS CAR T-CELL CONSTRUCT ANATOMY
         
  [ Extracellular Binding ]  ===> scFv (Single-chain variable fragment: VL - Linker - VH)
  [ Hinge & Transmembrane ] ===> CD8 / CD28 Hinge + Transmembrane Domain
  [ Costimulatory Domain ]  ===> 4-1BB (CD137) [Prolonged persistence] OR CD28 [Rapid, potent killing]
  [ T-Cell Activation ]     ===> CD3-zeta (ITAMs for intracellular signaling cascade)
```

### A. Autologous Anti-CD19 CAR T-Cell Therapies
* **Target & Mechanism**: Targets **CD19**, expressed from early B-cell development through mature B cells (lost on plasma cells). Infused T cells recognize CD19 on normal and malignant B cells, causing tumor lysis and profound **on-target, off-tumor B-cell aplasia** (managed with life-long monthly IVIG).

| Generic Name (INN) | Trade / Brand Name | Costimulatory Domain | Approved Clinical Indications | Medical Exam Pearls & Differentiators |
| :--- | :--- | :--- | :--- | :--- |
| **Tisagenlecleucel** | **Kymriah** | **4-1BB (CD137)** | • **Pediatric and Young Adult (up to 25 yr) Relapsed/Refractory B-cell ALL**<br>• Relapsed/refractory DLBCL | • **First CAR-T therapy approved by the FDA (2017)**.<br>• 4-1BB endodomain promotes mitochondrial biogenesis, oxidative phosphorylation, and long-term T-cell persistence/memory. |
| **Axicabtagene ciloleucel** | **Yescarta** | **CD28** | • **Large B-Cell Lymphoma (LBCL / DLBCL)** (second-line for early relapse/refractory, and $\ge 2$ lines)<br>• Relapsed Follicular Lymphoma | • CD28 endodomain promotes rapid aerobic glycolysis, leading to peak effector expansion and rapid cytoreduction, but faster T-cell exhaustion and higher peak CRS/ICANS rates. |
| **Brexucabtagene autoleucel** | **Tecartus** | **CD28** | • **Mantle Cell Lymphoma (MCL)** (relapsed/refractory)<br>• Adult Relapsed/Refractory **B-cell ALL** | • Incorporates a specialized T-cell enrichment step during manufacturing to remove circulating CD19+ circulating lymphoblasts, preventing premature activation and exhaustion. |
| **Lisocabtagene maraleucel** | **Breyanzi** | **4-1BB** | • **Large B-Cell Lymphoma (DLBCL)**<br>• Relapsed Follicular Lymphoma, Mantle Cell Lymphoma<br>• **Chronic Lymphocytic Leukemia (CLL / SLL)** | • Formulated as a **defined 1:1 ratio of CD4+ to CD8+ CAR T cells**, conferring highly predictable pharmacokinetic expansion and lower rates of Grade 3+ CRS/ICANS. |

---

### B. Autologous Anti-BCMA CAR T-Cell Therapies
* **Target & Mechanism**: Targets **B-Cell Maturation Antigen (BCMA / TNFRSF17)**, universally expressed on late-stage B cells, plasmablasts, normal plasma cells, and multiple myeloma cells.

| Generic Name (INN) | Trade / Brand Name | Molecular Architecture | Approved Indications | Clinical Pearls & Exam Points |
| :--- | :--- | :--- | :--- | :--- |
| **Idecabtagene vicleucel** | **Abecma** | Single scFv targeting BCMA + **4-1BB** domain | • **Relapsed/Refractory Multiple Myeloma** ($\ge 2$ prior lines, including IMiD, PI, anti-CD38) | • First CAR-T approved for multiple myeloma (KarMMa trial). Causes pan-hypogammaglobulinemia. |
| **Ciltacabtagene autoleucel** | **Carvykti** | **Bivalent** CAR with 2 single-domain camelid antibodies binding 2 distinct BCMA epitopes + **4-1BB** | • **Relapsed/Refractory Multiple Myeloma** (as early as 1 prior line with lenalidomide resistance) | • CARTITUDE trials: Unprecedented overall response rates ($>97\%$).<br>• **Unique Board Toxicity**: Delayed neurotoxicity presenting as **Parkinsonism / Movement Disorders** (cranial nerve palsies, tremors, micrographia) appearing weeks after initial CRS resolution. |

---

### C. Tumor-Infiltrating Lymphocytes (TIL) & Engineered TCR-T Therapies

| Generic Name (INN) | Trade / Brand Name | Modality & Molecular Target | Approved Clinical Indications | Critical Medical Exam Takeaways |
| :--- | :--- | :--- | :--- | :--- |
| **Lifileucel** | **Amtagvi** | **Autologous Tumor-Infiltrating Lymphocytes (TIL)** | • Unresectable or Metastatic **Melanoma** previously treated with anti-PD-1 $\pm$ BRAF/MEK inhibitors | • **First FDA-approved cellular therapy for solid tumors (2024)**.<br>• **Regimen Protocol**: Tumor resection $\to$ isolation/ex vivo expansion of polyclonal TILs $\to$ non-myeloablative lymphodepletion (Cy/Flu) $\to$ single TIL infusion $\to$ **up to 6 doses of high-dose systemic Aldesleukin (IL-2)** to support in vivo expansion. High toxicity burden! |
| **Afamitresgene autoleucel** (Afami-cel) | **Tecelra** | **Engineered TCR-T**: Autologous T cells transduced with TCR targeting **MAGE-A4** peptide restricted by **HLA-A\*02** | • Advanced, unresectable or metastatic **Synovial Sarcoma** (post-chemotherapy, HLA-A\*02-positive, MAGE-A4-positive) | • **First TCR-T therapy approved for solid tumors**.<br>• Unlike CAR-T (which only recognizes surface antigens), TCR-T recognizes intracellular proteins processed and presented on MHC Class I molecules. |

---

## 4. Costimulatory Receptor Agonists & Immune Modulators

While checkpoint inhibitors block negative signals ("cut the brakes"), costimulatory agonists stimulate activating pathways ("press the gas pedal") on antigen-presenting cells and effector T cells.

```
       ANTIGEN PRESENTING CELL                                    T CELL / NK CELL
      +------------------------+                               +--------------------+
      |  CD40                  | <==== [CD40 Agonist] ======== | CD40L (CD154)      |
      |  4-1BBL (CD137L)       | ======== [4-1BB Agonist] ===> | 4-1BB (CD137)      |
      |  OX40L (CD252)         | ======== [OX40 Agonist] ====> | OX40 (CD134)       |
      |  GITRL                 | ======== [GITR Agonist] ====> | GITR (CD357)       |
      +------------------------+                               +--------------------+
```

| Drug / Agent Name | Molecular Target & Mechanism | Development Status | Key Clinical / Exam Context |
| :--- | :--- | :--- | :--- |
| **Sotigalimab** (APX005M) | **CD40 Agonist**: Monoclonal antibody stimulating CD40 on dendritic cells and macrophages, licensing them to prime naive T cells and convert immunosuppressive "cold" tumors into "hot" inflamed tumors. | Clinical Trials (Phase I/II: Pancreatic cancer, melanoma) | Evaluated in the **Padron-iAtlas** cohort in metastatic pancreatic ductal adenocarcinoma in combination with nivolumab and chemotherapy. |
| **Urelumab** (BMS-663513) | **4-1BB (CD137) Agonist**: Fully human IgG4 activating 4-1BB costimulatory signaling on CD8+ T cells and NK cells. | Clinical Trials | **Classic Pharmacological Toxicity Lesson**: Severe, dose-limiting hepatotoxicity (autoimmune hepatitis) caused by liver macrophage cross-linking; modern 4-1BB agonists are engineered to require tumor-localized cross-linking. |
| **Utomilumab** (PF-05082566) | **4-1BB (CD137) Agonist**: Humanized IgG2 agonist; weaker clustering ability, markedly lower hepatotoxicity. | Clinical Trials (Solid tumors) | Better tolerated but lower single-agent clinical activity than urelumab. |
| **Imiquimod** (*Aldara*, *Zyclara*) | **Toll-Like Receptor 7 (TLR7) Agonist**: Small molecule imidazoquinoline stimulating TLR7 in plasmacytoid dendritic cells, inducing massive production of Interferon-alpha, IL-12, and TNF-alpha. | **FDA Approved** (Topical formulation) | • **Superficial Basal Cell Carcinoma (BCC)**<br>• Actinic Keratosis<br>• External Genital Warts (Condylomata acuminata). |
| **Vidutolimod** (CMP-001) | **Toll-Like Receptor 9 (TLR9) Agonist**: Unmethylated CpG-A oligodeoxynucleotide packaged inside a virus-like particle (bacteriophage Q-beta capsid). | Clinical Trials (Intratumoral injection in melanoma) | Injected intratumorally to stimulate plasmacytoid DCs, trigger systemic IFN-$\alpha$ release, and overcome anti-PD-1 resistance. |

---

## 5. Cytokines & Engineered Cytokine Receptor Agonists

Cytokines serve as critical autocrine and paracrine signaling molecules coordinating the magnitude and duration of anti-tumor immune responses.

| Generic Name (INN) | Trade / Brand Name | Mechanism of Action & Receptor Subunits | Approved Indications | Classic Medical Exam Toxicities & Pearls |
| :--- | :--- | :--- | :--- | :--- |
| **Aldesleukin** | **Proleukin** | **Recombinant Human Interleukin-2 (rhIL-2)**:<br>Binds the trimeric IL-2 receptor ($\alpha$ [CD25], $\beta$ [CD122], $\gamma$ [CD132]) on effector T cells, NK cells, and regulatory T cells. | • Metastatic **Renal Cell Carcinoma**<br>• Metastatic **Melanoma** (historic cure in a small subset [~7%] of patients) | • **Vascular / Capillary Leak Syndrome (CLS)**: IL-2 stimulates endothelial nitric oxide production and extravasation of fluid into interstitial spaces, causing **severe hypotension, diffuse non-cardiogenic pulmonary edema, dramatic weight gain, anuria, and shock**.<br>• Requires intensive care unit (ICU) monitoring during administration. |
| **Nogapendekin alfa inbakicept** | **Anktiva** (N-803) | **Interleukin-15 (IL-15) Superagonist**:<br>IL-15 mutant (N72D) bound to the dimeric IL-15 receptor alpha (IL-15R$\alpha$) sushi domain fused to human IgG1 Fc. Trans-presents IL-15 directly to the $\beta\gamma$ receptor on CD8+ T and NK cells **without activating immunosuppressive CD25+ Tregs**. | • **Non-Muscle Invasive Bladder Cancer (NMIBC)** with carcinoma in situ (CIS) that is BCG-unresponsive (FDA approved 2024, combined with intravesical BCG) | • **Exam Distinction vs. IL-2**: IL-15 does not bind CD25/IL-2R$\alpha$, avoiding Treg activation and significantly reducing capillary leak syndrome while driving potent CD8+ and NK cell memory expansion. |
| **Interferon alfa-2b**<br>*Peginterferon alfa-2b* | **Intron A**<br>*Sylatron* | **Recombinant Type I Interferon**:<br>Binds IFNAR1/IFNAR2, activating the JAK-STAT pathway (STAT1/2), upregulating MHC Class I expression, promoting NK/T-cell activation, and inhibiting tumor cell division. | • Adjuvant therapy for high-risk resected **Melanoma**<br>• Hairy Cell Leukemia<br>• AIDS-related Kaposi Sarcoma<br>• Follicular Lymphoma | • **High-Yield Side Effects**: Severe flu-like syndrome (fever, chills, myalgias, arthralgias), myelosuppression (neutropenia, thrombocytopenia), transaminitis, and **severe neuropsychiatric depression / suicide risk**. |
| **Sargramostim** | **Leukine** | **Recombinant GM-CSF** (Granulocyte-Macrophage Colony-Stimulating Factor): Stimulates proliferation and differentiation of granulocytes and macrophages; enhances dendritic cell antigen presentation. | • Accelerating myeloid reconstitution after autologous/allogeneic bone marrow transplantation<br>• Combined with Dinutuximab for neuroblastoma | • **Side Effects**: First-dose reaction, bone pain, fluid retention, dyspnea, and peripheral edema. |

---

## 6. Therapeutic Cancer Vaccines & Oncolytic Viruses

Cancer vaccines expose the immune system to tumor antigens along with strong adjuvants to induce de novo tumor-specific T-cell immunity. Oncolytic viruses selectively replicate within and lyse malignant cells, releasing damage-associated molecular patterns (DAMPs) and neoantigens in situ.

```
       ONCOLYTIC VIRUS (e.g., T-VEC) INJECTION
       
       [ Engineered Virus: HSV-1 (Delta ICP34.5 / ICP47 + GM-CSF) ]
                                |
                                v
               Infects Tumor Cell (Defective PKR / IFN pathway)
                                |
                   Replication & Lytic Destruction
                                |
          +---------------------+---------------------+
          |                                           |
          v                                           v
   Release of GM-CSF                     Release of Tumor Neoantigens
          |                                           |
          +---------------------+---------------------+
                                |
                                v
             Recruitment & Activation of Dendritic Cells
                                |
                                v
        Systemic Anti-Tumor CD8+ T Cell Priming (Abscopal Effect)
```

| Generic / Product Name | Trade Name | Modality & Engineering | Approved Indications | Medical Exam Board Highlights |
| :--- | :--- | :--- | :--- | :--- |
| **Talimogene laherparepvec** | **T-VEC** (*Imlygic*) | **Oncolytic Herpes Simplex Virus Type 1 (HSV-1)**:<br>• Deleted for **$\gamma 34.5$ (ICP34.5)** (prevents replication in normal cells with intact PKR antiviral pathways)<br>• Deleted for **ICP47** (restores MHC Class I antigen presentation)<br>• Inserted human **GM-CSF** gene. | • Local treatment of unresectable cutaneous, subcutaneous, and nodal lesions in patients with **Melanoma** recurrent after initial surgery. | • Injected **directly into accessible melanoma lesions**.<br>• Can induce regression of uninjected, distant visceral metastases through systemic immune priming (an **abscopal-like immune effect**).<br>• Healthcare workers must avoid accidental inoculation (risk of herpetic lesions). |
| **Nadofaragene firadenovec** | **Adstiladrin** | **Non-Replicating Adenoviral Vector (Ad5)**:<br>Delivers the human **Interferon alfa-2b (IFN$\alpha$-2b)** gene into urothelial bladder cells. Formulated with the excipient Syn3 to enhance viral transduction. | • High-risk, BCG-unresponsive **Non-Muscle Invasive Bladder Cancer (NMIBC)** with CIS $\pm$ high-grade Ta/T1 papillary disease. | • Instilled intravesically once every 3 months.<br>• Transduced urothelial cells turn into local "bio-factories" secreting high concentrations of IFN$\alpha$-2b, causing tumor cell death without systemic interferon toxicity. |
| **Sipuleucel-T** | **Provenge** | **Autologous Cellular Cancer Vaccine**:<br>Patient's peripheral blood mononuclear cells (including APCs) harvested via leukapheresis, incubated ex vivo with a fusion protein of **Prostatic Acid Phosphatase (PAP)** and **GM-CSF**, and reinfused over 3 cycles. | • Asymptomatic or minimally symptomatic **Metastatic Castration-Resistant Prostate Cancer (mCRPC)**. | • **High-Yield Clinical Trial Landmark (IMPACT Trial)**: Prolongs overall survival **without changing serum PSA levels or causing objective tumor shrinkage on CT/bone scans**. |
| **Bacillus Calmette-Guérin (BCG)** | *TICE BCG*, *TheraCys* | **Live Attenuated Mycobacterium bovis**:<br>Instilled directly into the bladder via catheter. Triggers a massive localized acute and chronic inflammatory reaction (macrophages, Th1 cells, IFN-$\gamma$, IL-2, IL-12), leading to immune destruction of superficial neoplastic urothelium. | • Primary intravesical therapy and prophylaxis for **High-Risk Non-Muscle Invasive Bladder Cancer (NMIBC)**. | • **Absolute Contraindications**: Gross hematuria (risk of systemic mycobacterial dissemination / "BCG sepsis"), traumatic catheterization, active tuberculosis, or severe immunosuppression.<br>• **Management of BCG Sepsis**: Immediate cessation + antitubercular therapy (Isoniazid, Rifampin, Ethambutol; note: *M. bovis* is intrinsically resistant to Pyrazinamide!). |
| **mRNA-4157 / V940**<br>*(Investigational)* | *(Individualized Neoantigen Vaccine)* | **Personalized Lipid-Nanoparticle mRNA Vaccine** encoding up to 34 patient-specific somatic neoepitopes identified by next-generation whole-exome sequencing. | • Adjuvant treatment of resected high-risk melanoma (Phase III, combined with Pembrolizumab). | Represents the frontier of precision oncology vaccines: trains T cells against private mutation antigens unique to an individual patient's tumor. |

---

## 7. Innate Checkpoints & ADCC-Mediating Monoclonals

Certain monoclonal antibodies target tumor-associated surface antigens and drive cytotoxicity primarily through host innate immunity via **Antibody-Dependent Cellular Cytotoxicity (ADCC)** mediated by CD16 (Fc$\gamma$RIIIa) on NK cells, **Antibody-Dependent Cellular Phagocytosis (ADCP)** mediated by macrophages, and **Complement-Dependent Cytotoxicity (CDC)**.

```
                    NK CELL (ADCC)                         MACROPHAGE (ADCP)
               +----------------------+                 +----------------------+
               | CD16a (Fc-gamma-RIIIa)|                |  Fc-gamma-RI / RIIa  |
               +----------------------+                 +----------------------+
                          |                                        |
                          v                                        v
                 [ Fc Region of mAb ]                    [ Fc Region of mAb ]
                          |                                        |
                          v                                        v
               [ Fab Binding Antigen ]                  [ Fab Binding Antigen ]
                          |                                        |
               +----------------------+                 +----------------------+
               |      Tumor Cell      |                 |      Tumor Cell      |
               | (CD20, CD38, or GD2) |                 | (e.g. CD47 Blockade) |
               +----------------------+                 +----------------------+
```

| Generic Name (INN) | Trade / Brand Name | Molecular Target & Mechanism | Approved Clinical Indications | Medical Exam Pearls & Side Effects |
| :--- | :--- | :--- | :--- | :--- |
| **Rituximab** | **Rituxan** | **Anti-CD20** chimeric IgG1:<br>Triggers profound ADCC, CDC, and apoptosis of CD20+ B cells. | • Non-Hodgkin Lymphomas (DLBCL, Follicular)<br>• Chronic Lymphocytic Leukemia (CLL)<br>• Autoimmune: Rheumatoid arthritis, GPA (Wegener's), Pemphigus vulgaris | • **Boxed Warnings**: Severe infusion reactions (premedicate with acetaminophen and diphenhydramine), Tumor Lysis Syndrome, Hepatitis B virus reactivation (mandatory screening before therapy!), and **Progressive Multifocal Leukoencephalopathy (PML)** caused by JC virus reactivation. |
| **Daratumumab** | **Darzalex** | **Anti-CD38** fully human IgG1:<br>Targets CD38 transmembrane glycoprotein on plasma cells. Drives ADCC, ADCP, CDC, and direct apoptosis. | • **Multiple Myeloma** (1st-line and relapsed/refractory; combined with IMiDs and proteasome inhibitors). | • **Classic Blood Bank Pitfall**: CD38 is weakly expressed on normal red blood cells. Daratumumab binds RBCs and causes **false-positive indirect Coombs tests (antibody screens)**, interfering with crossmatching. Blood typing and baseline Coombs testing must be performed *before* starting daratumumab. |
| **Isatuximab** | **Sarclisa** | **Anti-CD38** chimeric IgG1 binding a unique, distinct CD38 epitope. Direct ectoenzyme inhibition. | • Relapsed/Refractory Multiple Myeloma (+ pomalidomide/dex or carfilzomib/dex). | • Alternative CD38 inhibitor with minimal cross-resistance. Also interferes with RBC cross-matching. |
| **Elotuzumab** | **Empliciti** | **Anti-SLAMF7 (CS1)** humanized IgG1:<br>Binds SLAMF7 on myeloma cells (marking them for NK-cell mediated ADCC) AND binds SLAMF7 on NK cells (directly activating NK-cell cytotoxicity). | • Relapsed/Refractory Multiple Myeloma (+ lenalidomide/dex or pomalidomide/dex). | • Dual-action immunotherapy that directly boosts NK-cell intrinsic activation while simultaneously tagging myeloma cells. |
| **Dinutuximab** | **Unituxin** | **Anti-GD2** chimeric IgG1:<br>Binds the disialoganglioside GD2, heavily expressed on neuroblastoma cells and neuroectodermal tissues. | • Pediatric **High-Risk Neuroblastoma** (in maintenance with GM-CSF, IL-2, and isotretinoin). | • **Extreme Toxicity (Board Classic)**: GD2 is also expressed on peripheral nociceptive nerve fibers. Infusion induces **intense, excruciating neuropathic pain** requiring **continuous IV opioid infusions** before, during, and after therapy. |
| **Mogamulizumab** | **Poteligeo** | **Anti-CCR4** defucosylated humanized IgG1:<br>Selectively binds Chemokine Receptor 4 on T helper type 2 (Th2) and regulatory T cells (Tregs). Enhanced ADCC due to non-fucosylated Fc. | • **Mycosis Fungoides** and **Sézary Syndrome** (Cutaneous T-Cell Lymphomas [CTCL]) refractory to prior systemic therapies. | • Depletes immunosuppressive Tregs in the skin microenvironment. High risk of severe drug rashes and Stevens-Johnson syndrome (SJS). |
| **Magrolimab**<br>*(Investigational)* | *(Anti-CD47)* | **Anti-CD47 Blockade ("Don't Eat Me" Signal)**:<br>CD47 binds **SIRP$\alpha$** on macrophages to suppress phagocytosis. Magrolimab blocks CD47, unleashing macrophage-mediated phagocytosis. | • Investigated in Myelodysplastic Syndrome (MDS), Acute Myeloid Leukemia (AML), and solid tumors. | • **Key Toxicity**: Aging red blood cells rely on CD47 to prevent premature splenic clearance. Anti-CD47 antibodies frequently cause **hemagglutination, reticulocytosis, and acute hemolytic anemia**. |

---

## 8. Immunomodulatory Drugs (IMiDs) & Small-Molecule Checkpoints

Immunomodulatory drugs (IMiDs) are orally bioavailable small-molecule compounds that exert direct anti-proliferative, anti-angiogenic, and immunomodulatory effects via ubiquitin ligase reprogramming.

```
       [ IMiD: Thalidomide / Lenalidomide / Pomalidomide ]
                                |
                                v
               Binds Intracellular Receptor: CEREBLON (CRBN)
                                |
                                v
             Alters Substrate Specificity of CRL4-CRBN E3 Ubiquitin Ligase
                                |
          +---------------------+---------------------+
          |                                           |
          v                                           v
  Ubiquitination & Proteasomal                 Destabilization of Pro-Survival
  Degradation of IKZF1 (Ikaros) &              Transcription Factors (IRF4, MYC)
  IKZF3 (Aiolos) Zinc-Finger Proteins                         |
          |                                           Anti-Myeloma Cytotoxicity
          v
   Induces IL-2 Gene Transcription in T & NK Cells
          |
   Spurs Cytotoxic T-Cell & NK-Cell Immune Activation
```

| Generic Name (INN) | Trade / Brand Name | Primary Indications | High-Yield Medical Exam Board Points |
| :--- | :--- | :--- | :--- |
| **Thalidomide** | **Thalomid** | • Newly diagnosed Multiple Myeloma (+ dexamethasone)<br>• Erythema Nodosum Leprosum (ENL; acute cutaneous manifestations) | • **Historic Teratogen**: Causes severe **phocomelia** (flipper-like limbs) and amelia when taken during weeks 3–8 of pregnancy. Requires strict registration programs (STEPS / REMS).<br>• High incidence of **peripheral neuropathy**, constipation, and sedation. |
| **Lenalidomide** | **Revlimid** | • **Multiple Myeloma** (maintenance and combination regimens)<br>• **Myelodysplastic Syndrome (MDS)** associated with a **deletion 5q [del(5q)]** cytogenetic abnormality<br>• Mantle Cell Lymphoma, Follicular Lymphoma | • Significantly more potent than thalidomide with lower sedation and neuropathy, but significantly higher **myelosuppression (neutropenia, thrombocytopenia)**.<br>• **Black Box Warning**: Markedly increased risk of **Deep Vein Thrombosis (DVT) and Pulmonary Embolism (PE)** when combined with dexamethasone; mandatory thromboprophylaxis (aspirin or therapeutic LMWH) is indicated! |
| **Pomalidomide** | **Pomalyst** | • **Relapsed / Refractory Multiple Myeloma** (in patients who have received $\ge 2$ prior therapies including lenalidomide and a proteasome inhibitor)<br>• AIDS-related Kaposi Sarcoma | • Most potent approved IMiD; active even in lenalidomide-refractory disease. Also requires mandatory venous thromboembolism (VTE) prophylaxis. |

---

## 9. Master Medical Board Pharmacology & Toxicity Guide

Immune-related Adverse Events (irAEs), Cytokine Release Syndrome (CRS), and Neurotoxicity (ICANS) represent the most commonly tested clinical scenarios on medical examinations.

### A. Immune-Related Adverse Events (irAEs) — Comparative Pathology

```
     ORGAN SYSTEM                    ANTI-CTLA-4 (Ipilimumab)                ANTI-PD-1 / ANTI-PD-L1 (Pembro / Nivo)
     ================================================================================================================
     Endocrine System                HYPOPHYSITIS (~10-15%)                  THYROIDITIS (~5-10%)
                                     (Headache, hypopituitarism)             (Transient hyper -> permanent hypo)
     ----------------------------------------------------------------------------------------------------------------
     Gastrointestinal                SEVERE ENTEROCOLITIS (~30-40%)          MILD-TO-MODERATE COLITIS (~1-5%)
                                     (Watery diarrhea, perforation)          (Can present late)
     ----------------------------------------------------------------------------------------------------------------
     Pulmonary                       Rare (<1%)                              AUTOIMMUNE PNEUMONITIS (~3-5%)
                                                                             (Cough, dyspnea, "ground-glass" on CT)
     ----------------------------------------------------------------------------------------------------------------
     Dermatology                     Pruritic maculopapular rash             Vitiligo, Lichenoid rashes
                                                                             (Vitiligo predicts good tumor response!)
     ----------------------------------------------------------------------------------------------------------------
     Cardiovascular (Combined)       AUTOIMMUNE MYOCARDITIS (<1% incidence, but 30-50% fatality!)
                                     (Elevated troponin, conduction blocks; immediate high-dose methylprednisolone)
```

#### Detailed Organ Manifestations:
1. **Endocrine**:
   - **Hypophysitis (Pituitaritis)**: Strongly associated with **Ipilimumab (anti-CTLA-4)** due to ectopic CTLA-4 expression on anterior pituitary endocrine cells. Presents with dull frontal headache, fatigue, libido loss, secondary adrenal insufficiency (low ACTH, low cortisol), and central hypothyroidism (low TSH, low free T4). Pituitary enlargement on MRI with stalk thickening. Requires life-long hormone replacement (levothyroxine, hydrocortisone).
   - **Thyroid Disorders**: Strongly associated with **anti-PD-1/PD-L1**. Typically begins as painless destructive autoimmune thyroiditis with transient hyperthyroidism (suppressed TSH, elevated free T4), followed weeks later by permanent primary hypothyroidism (elevated TSH, low free T4). *High-yield tip*: Check thyroid function tests (TSH, free T4) every cycle; do not discontinue immunotherapy for simple hypothyroidism—initiate levothyroxine replacement!
   - **Type 1 Diabetes Mellitus**: Rare ($1\%$) rapid-onset autoimmune destruction of pancreatic beta cells presenting as acute, fulminant diabetic ketoacidosis (DKA) with low/undetectable C-peptide. Requires immediate insulin therapy; rarely recovers.
2. **Gastrointestinal**:
   - **Enterocolitis**: Severe cramping, watery non-bloody or bloody diarrhea, hematochezia. Risk of bowel perforation (most common cause of death from ipilimumab in early trials).
3. **Hepatic**:
   - **Immune-Mediated Hepatitis**: Asymptomatic acute elevations in serum transaminases (AST, ALT) $\pm$ total bilirubin. Histology shows pan-lobular hepatitis with CD8+ T-cell infiltration.
4. **Pulmonary**:
   - **Autoimmune Pneumonitis**: Dry cough, exertional dyspnea, fever, chest pain, and hypoxia. Chest CT shows bibasilar patchy ground-glass opacities or cryptogenic organizing pneumonia (COP) patterns.
5. **Dermatologic**:
   - **Vitiligo**: Autoimmune destruction of normal melanocytes sharing antigens with melanoma cells (e.g., MART-1, gp100). **Classic board question**: The emergence of vitiligo during checkpoint blockade in metastatic melanoma is a **favorable prognostic biomarker directly correlated with superior progression-free and overall survival!**

---

### B. Step-by-Step irAE Clinical Management Algorithm

| CTCAE Toxicity Grade | Clinical Severity Definition | Primary Management Action | Corticosteroid Regimen & Specialty Interventions |
| :--- | :--- | :--- | :--- |
| **Grade 1** | Mild symptoms; asymptomatic or clinical/diagnostic observations only; intervention not indicated. | **Continue Immunotherapy** (Exception: neurological, hematologic, or cardiac toxicities). | Close clinical observation and symptomatic relief (e.g., loperamide for mild diarrhea, topical hydrocortisone for mild rash). |
| **Grade 2** | Moderate symptoms; minimal, local, or non-invasive intervention indicated; limiting age-appropriate instrumental ADL. | **Temporarily Withhold Immunotherapy**. | • Initiate **Oral Prednisone $0.5\text{–}1.0\text{ mg/kg/day}$** (or equivalent).<br>• When symptoms improve to $\le \text{Grade } 1$, **slowly taper corticosteroids over at least 4 to 6 weeks**.<br>• Resume immunotherapy once steroid dose is $\le 10\text{ mg}$ prednisone equivalent. |
| **Grade 3** | Severe or medically significant but not immediately life-threatening; hospitalization or prolongation of existing hospitalization indicated. | **Permanently Discontinue Immunotherapy** (With rare exceptions for reversible endocrine toxicities controlled on replacement). | • Hospitalize immediately.<br>• Initiate **High-Dose IV Methylprednisolone $1\text{–}2\text{ mg/kg/day}$**.<br>• Once stabilized, transition to oral steroids and taper over at least 4–6 weeks.<br>• Provide gastric protection (PPI) and PCP prophylaxis (trimethoprim-sulfamethoxazole) during prolonged steroid tapers. |
| **Grade 4** | Life-threatening consequences; urgent intervention indicated. | **Permanently Discontinue Immunotherapy**. | Immediate ICU admission + IV methylprednisolone ($1\text{–}2\text{ mg/kg/day}$ or pulse doses of $1000\text{ mg/day}$ for myocarditis/encephalitis). |

#### Critical Refractory irAE Second-Line Immunosuppression (HIGH-YIELD BOARD EXAM TRAP!):
* **Steroid-Refractory Colitis** (no improvement after 48–72 hours of IV methylprednisolone):
  - Administer **Infliximab** (anti-TNF-$\alpha$ monoclonal antibody, $5\text{ mg/kg}$).
  - Alternative: **Vedolizumab** (gut-selective $\alpha_4\beta_7$ integrin antagonist; preferred in patients where systemic TNF inhibition is contraindicated).
* **Steroid-Refractory Hepatitis**:
  - **CRITICAL CONTRAINDICATION**: **DO NOT USE INFLIXIMAB!** Infliximab carries an FDA black box warning for drug-induced liver injury and can precipitate acute liver failure in patients with active hepatitis.
  - **CORRECT SECOND-LINE CHOICE**: Administer **Mycophenolate Mofetil (MMF)** ($500\text{–}1000\text{ mg}$ orally twice daily).

---

### C. Cytokine Release Syndrome (CRS) vs. Immune Effector Cell Neurotoxicity (ICANS)

CRS and ICANS are life-threatening systemic toxicities occurring after CAR T-cell infusions and BiTE therapies.

```
+---------------------------------------------------------------------------------------------------+
|                                 CRS vs. ICANS COMPARISON TABLE                                    |
+------------------------------+----------------------------------+---------------------------------+
| CLINICAL FEATURE             | CYTOKINE RELEASE SYNDROME (CRS)  | ICANS (NEUROTOXICITY)           |
+------------------------------+----------------------------------+---------------------------------+
| Peak Onset Time              | Early (Days 1 to 5 post-infusion)| Delayed (Days 4 to 10)          |
| Key Mediating Cytokines      | IL-6, IL-1, IFN-gamma, TNF-alpha | IL-1, IL-6 crossing BBB; GM-CSF |
| Cardinal Clinical Signs      | 1. High Fever (temperature >= 38C) | 1. Expressive Aphasia / Dysphasia|
|                              | 2. Hypotension (shock)           | 2. Handwriting Impairment (early)|
|                              | 3. Hypoxia                       | 3. Tremor, Confusion, Lethargy  |
|                              |                                  | 4. Seizures & Cerebral Edema    |
| First-Line Targeted Drug     | TOCILIZUMAB (IL-6R antagonist)   | DEXAMETHASONE (Corticosteroid)  |
| Role of Tocilizumab in ICANS | INEFFECTIVE / CONTRAINDICATED    | Does not penetrate BBB; blocks  |
|                              | AS MONOTHERAPY                   | peripheral IL-6R, causing serum |
|                              |                                  | IL-6 to spike and worsen CNS    |
|                              |                                  | toxicity!                       |
+------------------------------+----------------------------------+---------------------------------+
```

#### ASTCT Grading and Management Guide:
1. **CRS Grading**:
   - **Grade 1**: Fever ($\ge 38.0^\circ\text{C}$), NO hypotension, NO hypoxia. $\to$ Supportive care, antipyretics.
   - **Grade 2**: Fever + Hypotension responding to IV fluids OR Hypoxia requiring low-flow nasal cannula ($\le 40\%\text{ FiO}_2$). $\to$ **Tocilizumab** IV ($8\text{ mg/kg}$, max $800\text{ mg}$; may repeat every 8 hours, up to 4 doses) $\pm$ Dexamethasone.
   - **Grade 3**: Fever + Hypotension requiring a single high-dose vasopressor OR Hypoxia requiring high-flow nasal cannula ($>40\%\text{ FiO}_2$). $\to$ **Tocilizumab + High-dose Dexamethasone** ($10\text{–}20\text{ mg}$ IV every 6 hours).
   - **Grade 4**: Fever + Hypotension requiring multiple vasopressors OR Hypoxia requiring mechanical ventilation (CPAP, BiPAP, intubation). $\to$ **Tocilizumab + IV Methylprednisolone $1000\text{ mg/day}$** $\pm$ Siltuximab (anti-IL-6) or Anakinra (IL-1 receptor antagonist).
2. **ICANS Management**:
   - Monitored using the **10-point ICE (Immune Effector Cell-Associated Encephalopathy) score** (evaluates orientation, naming, following commands, writing a complete standard sentence, and counting backward from 100).
   - Deterioration in handwriting is the **earliest clinical warning sign** of impending neurotoxicity.
   - First-line therapy for Grade 2+ ICANS is **IV Dexamethasone** ($10\text{ mg}$ every 6 hours) because dexamethasone readily penetrates the blood-brain barrier.
   - Prophylactic **Levetiracetam (Keppra)** is routinely administered for seizure prophylaxis during cell therapy.

---

### D. Tumor Microenvironment Biomarkers & Predictive Signatures

```
     BIOMARKER                    TESTING METHOD                     CLINICAL INTERPRETATION / PREDICTIVE VALUE
     ===========================================================================================================
     MSI-H / dMMR                 IHC (MLH1, MSH2, MSH6, PMS2)       Tumor-Agnostic approval for Pembrolizumab /
                                  or PCR / NGS for microsatellites   Dostarlimab. Hundreds of frameshift mutations
                                                                     create foreign neoantigens.
     -----------------------------------------------------------------------------------------------------------
     Tumor Mutational Burden      Targeted NGS Panel / WES           >= 10 mutations / megabase (TMB-High) predicts
     (TMB-High)                                                      higher probability of durable ICI response.
     -----------------------------------------------------------------------------------------------------------
     PD-L1 Expression             Immunohistochemistry (IHC)         • TPS (Tumor Proportion Score): % viable tumor
                                  (Assays: 22C3, 28-8, SP142, SP263)   cells with membrane PD-L1 staining (used in NSCLC).
                                                                     • CPS (Combined Positive Score):
                                                                       [(PD-L1+ Tumor + Lymphocytes + Macrophages) /
                                                                       Total Tumor Cells] x 100 (used in Gastric,
                                                                       Cervical, HNSCC, TNBC).
     -----------------------------------------------------------------------------------------------------------
     Genomic Resistance           Whole-Exome Sequencing / RNA-seq   • Loss-of-function mutations in B2M or HLA-A/B/C
     Mutations                                                         abolish antigen presentation.
                                                                     • Inactivating mutations in JAK1 / JAK2 prevent
                                                                       IFN-gamma signaling and stop PD-L1 upregulation.
                                                                     • STK11 / KEAP1 co-mutations in KRAS-mutant NSCLC
                                                                       confer primary ICI resistance ("cold" microenvironment).
```

---

## 10. Master Rapid-Review Medical Exam Summary Table

| Medication (Generic) | Brand Name | Primary Target | Therapeutic Modality | Approved Tumors / Indications | Classic Board Exam Hallmarks |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Pembrolizumab** | Keytruda | PD-1 | Humanized mAb | Melanoma, NSCLC, MSI-H/dMMR, TNBC, RCC | First tumor-agnostic FDA approval based on genomic defect (MSI-H). |
| **Nivolumab** | Opdivo | PD-1 | Fully human mAb | Melanoma, NSCLC, RCC, cHL, Bladder, CRC | CheckMate-067 (+ Ipilimumab) set melanoma standard; IgG4 framework. |
| **Cemiplimab** | Libtayo | PD-1 | Fully human mAb | Cutaneous Squamous Cell Carcinoma (CSCC) | Standard of care for non-resectable locally advanced or metastatic CSCC. |
| **Dostarlimab** | Jemperli | PD-1 | Humanized mAb | dMMR Endometrial & Rectal Cancer | 100% complete clinical response rate in dMMR locally advanced rectal cancer. |
| **Atezolizumab** | Tecentriq | PD-L1 | Engineered humanized mAb | SCLC, Hepatocellular Carcinoma, NSCLC | Aglycosylated Fc (N298A); standard 1st-line in HCC combined with Bevacizumab. |
| **Durvalumab** | Imfinzi | PD-L1 | Fully human mAb | Stage III NSCLC, SCLC, Biliary Tract | Consolidation standard following definitive chemoradiation (PACIFIC trial). |
| **Avelumab** | Bavencio | PD-L1 | Fully human mAb | Urothelial maintenance post-platinum, MCC | Retains active IgG1 Fc domain capable of mediating direct NK-cell ADCC. |
| **Ipilimumab** | Yervoy | CTLA-4 | Fully human IgG1 | Melanoma, RCC, NSCLC, Mesothelioma, HCC | Priming phase inhibitor; high risk of autoimmune hypophysitis and enterocolitis. |
| **Tremelimumab** | Imjudo | CTLA-4 | Fully human IgG2 | Hepatocellular Carcinoma, NSCLC | Used as single priming dose in the STRIDE regimen with regular Durvalumab. |
| **Relatlimab** | Opdualag (+Nivo)| LAG-3 | Humanized mAb | Metastatic Melanoma | First-in-class LAG-3 checkpoint inhibitor; fixed-dose co-infusion with Nivolumab. |
| **Blinatumomab** | Blincyto | CD19 $\times$ CD3 | Bispecific T-cell Engager (BiTE) | B-cell Acute Lymphoblastic Leukemia (MRD+) | Short 2-hour half-life requires continuous 28-day IV infusion; black box for CRS. |
| **Tebentafusp** | Kimmtrak | gp100 $\times$ CD3 | TCR-anti-CD3 bispecific (ImmTAC) | Metastatic Uveal Melanoma | Strictly restricted to HLA-A*02:01-positive patients; first approved TCR therapy. |
| **Tarlatamab** | Imdelltra | DLL3 $\times$ CD3 | Bispecific T-cell Engager | Extensive-Stage Small Cell Lung Cancer | Targets Delta-like ligand 3 on neuroendocrine tumor cells. |
| **Teclistamab** | Tecvayli | BCMA $\times$ CD3 | Bispecific T-cell Engager | Relapsed/Refractory Multiple Myeloma | BCMA-directed T-cell redirector; requires IVIG for deep hypogammaglobulinemia. |
| **Talquetamab** | Talvey | GPRC5D $\times$ CD3 | Bispecific T-cell Engager | Relapsed/Refractory Multiple Myeloma | Targets GPRC5D; causes severe dysgeusia (loss of taste), skin peeling, nail dystrophy. |
| **Tisagenlecleucel**| Kymriah | CD19 | Autologous CAR-T (4-1BB) | Pediatric B-ALL, DLBCL | First CAR-T approved (2017); 4-1BB domain promotes oxidative persistence. |
| **Axicabtagene** | Yescarta | CD19 | Autologous CAR-T (CD28) | Large B-Cell Lymphoma (DLBCL), FL | CD28 domain promotes rapid effector expansion but higher peak CRS/ICANS rates. |
| **Ciltacabtagene** | Carvykti | BCMA | Bivalent CAR-T (4-1BB) | Multiple Myeloma | 2 camelid nanobodies; landmark response rates; delayed Parkinsonian neurotoxicity. |
| **Lifileucel** | Amtagvi | Polyclonal TIL | Autologous Tumor-Infiltrating T Cells | Advanced Metastatic Melanoma | First cell therapy approved for solid tumors; requires high-dose IL-2 post-infusion. |
| **Afami-cel** | Tecelra | MAGE-A4 $\times$ HLA-A*02 | Engineered TCR-T Cell Therapy | Advanced Synovial Sarcoma | Transgenic TCR recognizing intracellular MAGE-A4 peptide on MHC-I. |
| **Aldesleukin** | Proleukin | IL-2 | Recombinant Cytokine | Metastatic RCC, Melanoma | High-dose IL-2 causes severe capillary leak syndrome (hypotension, anasarca, shock). |
| **Nogapendekin** | Anktiva | IL-15 superagonist | Engineered Cytokine Complex | BCG-Unresponsive Bladder CIS | Selectively stimulates CD8/NK cells without activating CD25+ immunosuppressive Tregs. |
| **T-VEC** | Imlygic | HSV-1 + GM-CSF | Oncolytic Virus | Recurrent / Inoperable Melanoma | Injected directly into skin lesions; can stimulate systemic abscopal antitumor immunity. |
| **Nadofaragene** | Adstiladrin | Ad5-IFN$\alpha$-2b | Non-replicating Viral Vector | BCG-Unresponsive Bladder CIS | Instilled intravesically every 3 months; creates local urothelial interferon factories. |
| **Sipuleucel-T** | Provenge | PAP-GM-CSF | Autologous Cellular Vaccine | Metastatic Castration-Resistant Prostate Ca | Extends overall survival without changing serum PSA or objective tumor dimensions. |
| **BCG** | TICE BCG | *M. bovis* | Live Attenuated Mycobacteria | Non-Muscle Invasive Bladder Cancer | Live intravesical immunotherapy; contraindicated in gross hematuria (sepsis risk). |
| **Rituximab** | Rituxan | CD20 | Chimeric mAb | DLBCL, Follicular Lymphoma, CLL | Drives ADCC/CDC; boxed warnings for HBV reactivation, infusion shock, and PML. |
| **Daratumumab** | Darzalex | CD38 | Fully human IgG1 mAb | Multiple Myeloma | Binds CD38 on RBCs causing false-positive Coombs tests and crossmatch interference. |
| **Dinutuximab** | Unituxin | GD2 | Chimeric mAb | High-Risk Pediatric Neuroblastoma | Binds GD2 on peripheral nerves causing severe neuropathic pain (requires continuous IV opioids). |
| **Thalidomide** | Thalomid | Cereblon (CRBN) | Small Molecule IMiD | Multiple Myeloma, Leprosy (ENL) | Classic teratogen causing phocomelia; binds cereblon to degrade Ikaros/Aiolos. |
| **Lenalidomide** | Revlimid | Cereblon (CRBN) | Small Molecule IMiD | Multiple Myeloma, 5q- MDS | Higher myelosuppression and DVT/PE risk than thalidomide; mandatory thromboprophylaxis. |
| **Tocilizumab** | Actemra | IL-6 Receptor | Humanized mAb | Cytokine Release Syndrome (CRS) | First-line targeted therapy for CRS; ineffective as monotherapy for ICANS. |
