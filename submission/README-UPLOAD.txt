ENAQT paper v6 - upload guide (October 27, 2026)
=================================================

WHAT IS HERE
  enaqt_v6_arxiv_source.tar.gz   arXiv source: main.tex + figures/ (4 PNGs).
                                 The bibliography is inside main.tex, so no
                                 bibtex run is needed.
  main.pdf                       The compiled paper (pdfLaTeX, revtex4-2).
                                 Read this first.
  main.tex (repo root)           Same manuscript as in the tarball.
  abstract.txt                   Plain-text abstract for the arXiv form
                                 (1,871 characters; the limit is 1,920).

WHAT CHANGED SINCE v5
  - New central result: the optimal dephasing rate equals the
    nearest-neighbour energy mismatch, gamma* = eps_tot/(N-1), to within 2%.
    The N^-1.20 power law is now only a fit to this.
  - New exact result: the yield does not change when the bias is reversed,
    and the flat chain beats every funnel.
  - Robustness checks (other trap and loss rates, other total biases),
    bootstrap confidence intervals, matched-disorder comparison.
  - FMO at biological trap and loss rates: the enhancement falls to 1.14.
  - Corrected numbers (Table I, Table IV, factor 10.7 at N=15, BChl 7 ratio,
    optimizer profiles, localization lengths).
  - Prior work now credited, mainly Jacob et al. 2024 for the mismatch rule;
    19 new references; bibliography in citation order.

BEFORE YOU POST - CHECK THESE
  1. Author line: "A. Smith", affiliation "Ember Professional Research
     Division", email grandmasta1024@gmail.com (lines 39-41 of main.tex).
  2. FMO bath relaxation time. The paper uses 1/gamma_c = 50 fs for
     Ishizaki and Fleming, PNAS 106, 17255 (2009), the value most papers
     quote. One paper attributes 100 fs to the same source. If you can open
     the PNAS paper, check which value it uses. The text already says that
     100 fs would double the estimate, so no other change is needed.
  3. The GitHub repository's About line still mentions HEOM validation, a
     ~2.1N law and "disorder amplifies ENAQT". All are withdrawn. Change it
     on the repository page (gear icon next to About) before posting.

GETTING A PDF YOURSELF (optional)
  - Overleaf: New Project > Upload Project > choose the .tar.gz. Set the
    compiler to pdfLaTeX and click Recompile.
  - MiKTeX on Windows: unpack the tarball, open a terminal in that folder,
    run "pdflatex main" twice. MiKTeX installs revtex4-2 on first use.

POSTING TO arXiv
  1. Sign in at arxiv.org (or register).
  2. First quant-ph submission: arXiv will ask for an endorsement. It gives
     you an endorsement code; send that code to the person who agreed to
     endorse you.
  3. Start a new submission:
     - License: CC BY 4.0 is the usual choice.
     - Primary category: quant-ph.
     - Optional cross-list: physics.bio-ph.
  4. Upload enaqt_v6_arxiv_source.tar.gz as is. Check the PDF preview.
  5. Metadata:
     - Title: Environment-Assisted Quantum Transport in Open Quantum Chains:
       Optimal Dephasing, Disorder Statistics, and the Limits of the
       Enhancement Ratio
     - Abstract: paste abstract.txt.
     - Comments: "22 pages, 4 figures. Code and data:
       https://github.com/GrandMastaShake/enaqt-open-quantum-chains".
  6. Submit. It appears after the next daily announcement.

ABOUT bioRxiv
  bioRxiv does not accept papers already posted to arXiv. It wants a single
  PDF, and this paper is mainly physics. So use arXiv, or bioRxiv instead of
  arXiv, not both.

AFTER POSTING
  Merge the v6 pull request on GitHub so the repo matches the paper.
