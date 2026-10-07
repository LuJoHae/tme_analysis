# Cancer Immunotherapy — Medical Exam Anki Flashcard Deck

This folder contains the complete, high-yield spaced repetition Anki flashcard deck derived directly from the [`docs/medical_exam_immunotherapy_study_guide.md`](../medical_exam_immunotherapy_study_guide.md).

---

## Available Formats

1. **Native Anki Package (`immunotherapy_medical_exam.apkg`)**:
   - Double-clickable package file.
   - Pre-configured with the deck name **`Cancer_Immunotherapy::Medical_Board_Review`**.
   - Built-in card styling (color-coded badges for drugs, toxicities, targets, and exam pearls).
   - Works immediately on Anki Desktop (macOS, Windows, Linux).

2. **Universal Text Import Deck (`immunotherapy_medical_exam_deck.txt`)**:
   - Standard tab-delimited Anki import format with `#separator:tab`, `#html:true`, and `#tags column:3`.
   - Compatible with Anki Desktop, AnkiMobile (iOS), AnkiDroid (Android), and AnkiWeb.

---

## How to Import

### Option 1: Native `.apkg` File (Recommended — Fastest)
1. Open the **Anki** desktop application.
2. Double-click [`immunotherapy_medical_exam.apkg`](./immunotherapy_medical_exam.apkg) (or go to `File > Import...` and select `immunotherapy_medical_exam.apkg`).
3. Anki will automatically create the deck **`Cancer_Immunotherapy::Medical_Board_Review`** with all 73 styled cards and hierarchical tags.

### Option 2: Text Import File (`.txt`)
1. Open **Anki**.
2. Click **File > Import...** (or press `Ctrl+I` / `Cmd+I`).
3. Select [`immunotherapy_medical_exam_deck.txt`](./immunotherapy_medical_exam_deck.txt).
4. In the Import Dialog:
   - **Type**: `Basic` (or create a custom note type).
   - **Deck**: Select or create `Cancer_Immunotherapy::Medical_Board_Review`.
   - **Field 1**: Map to `Front`.
   - **Field 2**: Map to `Back`.
   - Ensure the checkbox **"Allow HTML in fields"** is checked.
5. Click **Import**.

---

## Deck Structure & Hierarchical Tags

All cards are tagged hierarchically to allow targeted cramming using Anki's **Custom Study / Filtered Decks**:

| Tag Hierarchy | Card Topics |
| :--- | :--- |
| `Immunotherapy::Pharmacology_Suffixes` | Monoclonal antibody suffixes (`-omab`, `-ximab`, `-zumab`, `-umab`, `-cel`, `-vec`, `-cept`) & IgG4 vs IgG1 Fc engineering. |
| `Immunotherapy::ICIs::PD1_PDL1` | Anti-PD-1 (Pembro, Nivo, Cemiplimab, Dostarlimab, Tislelizumab, Toripalimab, Retifanlimab) & Anti-PD-L1 (Atezolizumab, Durvalumab, Avelumab, Envafolimab). |
| `Immunotherapy::ICIs::CTLA4_LAG3` | Anti-CTLA-4 (Ipilimumab, Tremelimumab), LAG-3 (Relatlimab / Opdualag), TIGIT, Cadonilimab (PD-1/CTLA-4), Ivonescimab (PD-1/VEGF). |
| `Immunotherapy::BiTEs` | Bispecific T-cell engagers: Blinatumomab (CD19), Tebentafusp (gp100 HLA-A*02), Tarlatamab (DLL3), Teclistamab/Elranatamab (BCMA), Talquetamab (GPRC5D), Glofitamab/Mosunetuzumab/Epcoritamab (CD20). |
| `Immunotherapy::Cellular::CAR_T` | Anti-CD19 CAR-T (Kymriah, Yescarta, Tecartus, Breyanzi), Anti-BCMA CAR-T (Abecma, Carvykti), CD28 vs 4-1BB endodomains, lymphodepletion rationale. |
| `Immunotherapy::Cellular::TCR_TIL` | Lifileucel (TILs in melanoma + IL-2), Afami-cel (MAGE-A4 TCR-T in synovial sarcoma). |
| `Immunotherapy::Cytokines_Agonists` | Aldesleukin (IL-2 & capillary leak syndrome), Nogapendekin alfa (IL-15 superagonist), IFNα-2b, Imiquimod (TLR7), STING agonists, Nemvaleukin alfa. |
| `Immunotherapy::Vaccines_Oncolytic` | T-VEC (HSV-1 + GM-CSF abscopal effect), Adstiladrin, Sipuleucel-T (OS benefit without PSA change), intravesical BCG contraindications & BCG sepsis. |
| `Immunotherapy::Innate_IMiDs` | Rituximab (HBV/PML warnings), Daratumumab (Coombs RBC crossmatch interference), Dinutuximab (pain & opioids), CD47 "don't eat me", IMiDs (Cereblon, Thalidomide, Lenalidomide, Pomalidomide). |
| `Immunotherapy::Toxicities::irAEs` | Hypophysitis vs Thyroiditis, Vitiligo favorable prognosis, autoimmune myocarditis, corticosteroid tapering rules, Infliximab vs Mycophenolate Mofetil. |
| `Immunotherapy::Toxicities::CRS_ICANS` | Peak onset, pathophysiology, Tocilizumab for CRS vs Dexamethasone for ICANS, ICE scoring and handwriting loss. |
| `Immunotherapy::Biomarkers_Pharmacology` | MSI-H/dMMR & Lynch syndrome, TMB-H, TPS vs CPS, B2M/JAK1/JAK2/STK11 resistance mutations. |
| `Board_High_Yield` | Filter tag highlighting the most frequently tested board exam questions and pitfalls across all categories. |

---

## Recommended Study Strategy for Medical Exams

1. **Daily Spaced Repetition**: Study 20 new cards per day with unlimited reviews.
2. **Targeted Filtered Decks before Rounds / Exams**:
   - To review toxicities before clinical oncology rotations or emergency medicine wards:
     ```
     tag:Immunotherapy::Toxicities::irAEs or tag:Immunotherapy::Toxicities::CRS_ICANS
     ```
   - To review high-yield board favorites in the final 48 hours before the exam:
     ```
     tag:Board_High_Yield
     ```
