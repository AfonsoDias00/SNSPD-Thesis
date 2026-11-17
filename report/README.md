# LaTeX × TU Delft - Report/Thesis Template

This template aims to simplify and improve the (Xe)LaTeX report/thesis template by Delft University of Technology with the following three main design principles:

* **Simplicity First:** A class file that has been reduced by nearly 70% to simplify customization;
* **Effortless:** A careful selection of common packages to get started immediately;
* **Complete:** Ready-to-go when it comes to the document and file structure.

This template works with _pdfLaTeX_, _XeLaTeX_ and _LuaLaTeX_. In order to adhere to the TU Delft house style, either _XeLaTeX_ or _LuaLaTeX_ is required, as it supports TrueType and OpenType fonts. _BibLaTeX_ is used for the bibliography with as backend _biber_. Please visit https://dzwaneveld.github.io/report/ for the full documentation.

Cover | Title | Chapter
--- | --- | ---
<img src="https://dzwaneveld.github.io/images/report-template.jpg"> | <img src="https://dzwaneveld.github.io/images/report-template-title.jpg"> | <img src="https://dzwaneveld.github.io/images/report-template-chapter.jpg">

## Building this thesis

Requirements: a TeX distribution with `latexmk`, XeLaTeX and biber (e.g. TeX Live), and optionally VS Code with the *LaTeX Workshop* extension.

* **VS Code:** open this folder and save any `.tex` file (or run *LaTeX Workshop: Build LaTeX project*). The recipe is in `.vscode/settings.json`.
* **Command line:** run `latexmk report.tex` in this folder.

All build settings are in `.latexmkrc`. It compiles with XeLaTeX and writes the auxiliary files to `%TEMP%/latex-build/Thesis_report` (on Linux/WSL: `/tmp/latex-build/Thesis_report`), so this folder only receives `report.pdf` (and `report.synctex.gz` for jumping between source and PDF). If the build gets stuck, delete that temp folder.

**Fonts.** Body text is Libertinus Serif with matching maths (`newtxmath`), titles are Roboto Slab (TU Delft style); see `tudelft-report.cls`. Both come with TeX Live / the system, nothing to install.

## Layout of this thesis

* `mainmatter/01-introduction.tex` ... `07-conclusion.tex`: one file per chapter, included from `report.tex`.
* `appendix/a-validation.tex`, `b-parameters.tex`; `typesetting_examples.tex` is the template's cheat sheet for the macros (not built; uncomment its `\input` to see it).
* `figures/<chapter>/`: figures per chapter, copied from `../simulations/results/`. `\graphicspath` is `figures/`, so write `\includegraphics{04-results/name}`.
* `report.bib`: entries marked `% VERIFY` still need their details checked.
* `MScAP_MTP_Preparation_Report__Afonso.pdf`: the preparation report (6 Oct 2026), for reference.

## Documentation (Abridged)

As a report/thesis is generally a substantial document, the chapters and appendices have been separated into different files and folders for convenience. The folders are based on the three parts in the document: the frontmatter, mainmatter and appendix. All files are inserted in the main file, `report.tex`, using the `\input{filename}` command. The document class, which can be found in `tudelft-report.cls`, is based on the `book` class.

The template will automatically generate a cover when the `\makecover` command is used. The title, subtitle and author will also be present on the title page. To give greater flexibility over the title page, the layout is specified in `title-report.tex`. A title page for theses is also available: `title-thesis.tex`. Change the corresponding `\input{...}` command in the main file to switch. 

The bibliography has been set up in `report.tex` to allow for easy customization. It is included in the table of contents and renamed to 'References' using the `heading=bibintoc` and `title=References` options of the `\printbibliography` command respectively. If you would like to use a different `.bib` file, change the command `\addbibresource{report.bib}` accordingly. 

*→ Visit https://dzwaneveld.github.io/report/ for the full documentation.*

## License

This [report/thesis template](https://github.com/dzwaneveld/TU-Delft-Unofficial-Report-Template) by Daan Zwaneveld is licensed under [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/). No attribution is required in PDF outputs created using this template.
