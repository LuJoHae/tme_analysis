// Master entry point for the manuscript

#import "styles.typ": article_setup, title_block
#import "metadata.typ": paper_title, paper_authors, paper_affiliations, paper_abstract, paper_keywords

#show: article_setup

#title_block(
  paper_title,
  paper_authors,
  paper_affiliations,
  paper_abstract,
  paper_keywords
)

#columns(2, gutter: 1.5em)[
  #include "chapters/01_introduction.typ"
  #include "chapters/02_methods.typ"
  #include "chapters/03_results_single_cell.typ"
  #include "chapters/04_results_deconvolution.typ"
  #include "chapters/05_results_mutations.typ"
  #include "chapters/06_results_subclonal_vaf.typ"
  #include "chapters/07_results_tcell_myeloid.typ"
  #include "chapters/08_results_milo_vs_deconv.typ"
  #include "chapters/09_discussion.typ"

  #v(1.5em)
  #bibliography("references.bib", style: "nature")
]

#pagebreak()
#include "chapters/10_extended_data.typ"
