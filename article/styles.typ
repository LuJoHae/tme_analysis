// Centralized styling rules and template configuration for the article

#let article_setup(body) = {
  set page(
    paper: "a4",
    margin: (x: 1.5cm, y: 2.2cm),
    header: align(right)[#text(size: 8pt, fill: gray)[TME Deconvolution & Genomic Prediction]],
    footer: context [
      #align(center)[#text(size: 8pt)[#counter(page).display("1")]]
    ]
  )

  set text(
    font: ("Times New Roman", "Charter", "Palatino", "Georgia"),
    size: 9.5pt,
    fill: rgb("#2C3E50")
  )

  set par(justify: true, leading: 0.6em)

  show heading: it => [
    #set text(fill: rgb("#1F3A60"))
    #v(1.2em)
    #it
    #v(0.6em)
  ]

  show heading.where(level: 1): it => {
    block(width: 100%, below: 0.8em)[
      #set text(size: 11pt, weight: "bold")
      #it.body
    ]
  }

  show heading.where(level: 2): it => {
    block(width: 100%, below: 0.6em)[
      #set text(size: 10pt, weight: "bold", style: "italic")
      #it.body
    ]
  }

  show heading.where(level: 3): it => {
    block(width: 100%, below: 0.4em)[
      #set text(size: 9.5pt, weight: "bold")
      #it.body
    ]
  }

  show figure.caption: it => [
    #set text(size: 8pt, style: "italic")
    #it
  ]

  body
}

#let title_block(title, authors, affiliations, abstract, keywords) = {
  align(center)[
    #text(size: 15pt, weight: "bold", fill: rgb("#1F3A60"))[#title]
    #v(0.8em)
    #text(size: 11pt, style: "italic")[#authors]
    #v(0.4em)
    #text(size: 8.5pt, fill: gray)[#affiliations]
  ]

  v(1em)

  block(
    fill: rgb("#F8F9FA"),
    inset: 1.2em,
    radius: 4pt,
    stroke: 0.5pt + rgb("#E2E8F0"),
    width: 100%
  )[
    #align(center)[#text(weight: "bold", size: 10.5pt)[Abstract]]
    #v(0.4em)
    #text(size: 8.5pt)[#abstract]
    #if keywords != none [
      #v(0.5em)
      #text(size: 8pt, weight: "bold")[Keywords: ]
      #text(size: 8pt, style: "italic")[#keywords.join(", ")]
    ]
  ]

  v(0.8em)
}
