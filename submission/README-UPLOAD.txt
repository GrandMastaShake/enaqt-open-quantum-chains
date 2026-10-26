ENAQT paper v5 - upload guide (October 26, 2026)
=================================================

WHAT IS HERE
  enaqt_v5_arxiv_source.tar.gz   arXiv source: main.tex + figures/ (4 PNGs).
                                 The bibliography is inside main.tex, so no
                                 bibtex run is needed.
  main.tex                       Same manuscript as in the tarball.
  abstract.txt                   Plain-text abstract for the arXiv form.

  Everything is also on GitHub, branch paper-v5-submission-ready (draft PR #1):
  https://github.com/GrandMastaShake/enaqt-open-quantum-chains/pull/1

BEFORE YOU POST - CHECK THESE
  1. Author line: "A. Smith", affiliation "Ember Professional Research
     Division", email grandmasta1024@gmail.com. Edit lines 34-36 of main.tex
     if any of this should change.
  2. v5 changes the paper's main message. The old claims that HEOM validates
     the model, that disorder amplifies transport, and that the funnel is the
     optimal disorder are withdrawn, because the recomputed numbers do not
     support them. The new message is that the enhancement ratio is dominated
     by the coherent baseline. Read the abstract and Sec. III once.
  3. No PDF was built here (no LaTeX on the build machine). The source was
     checked by hand: braces, environments, labels, citations and figure
     paths all match. Look at the first PDF before announcing anything.

GETTING A PDF (pick one)
  - Overleaf: New Project > Upload Project > choose the .tar.gz. Set the
    compiler to pdfLaTeX and click Recompile.
  - MiKTeX on Windows: unpack the tarball, open a terminal in that folder,
    run "pdflatex main" twice. MiKTeX installs revtex4-2 on first use.
  - arXiv: the preview it shows after upload is the real build.

POSTING TO arXiv
  1. Sign in at arxiv.org (or register).
  2. First quant-ph submission: arXiv will ask for an endorsement. Ask
     someone who has posted to quant-ph to endorse you with the code arXiv
     gives you.
  3. Start a new submission:
     - License: CC BY 4.0 is the usual choice.
     - Primary category: quant-ph.
     - Optional cross-list: physics.bio-ph.
  4. Upload enaqt_v5_arxiv_source.tar.gz as is. Check the PDF preview.
  5. Metadata:
     - Title: Environment-Assisted Quantum Transport in Open Quantum Chains:
       Length Scaling, Disorder Statistics, and the Limits of the
       Enhancement Ratio
     - Abstract: paste abstract.txt.
     - Comments: "4 figures. Code and data:
       https://github.com/GrandMastaShake/enaqt-open-quantum-chains".
  6. Submit. It appears after the next daily announcement.

ABOUT bioRxiv
  bioRxiv does not accept papers already posted to arXiv. It wants a single
  PDF, and this paper is mainly physics. So use arXiv, or bioRxiv instead of
  arXiv, not both.

AFTER POSTING
  Merge PR #1 on GitHub so the repo matches the paper. CI only runs once
  the GitHub billing lock on your account is cleared.
