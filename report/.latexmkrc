# latexmk configuration for this thesis (used by VS Code LaTeX Workshop and the command line)

# Compile with XeLaTeX (required for the TU Delft house style fonts); biber runs automatically
$pdf_mode = 5;

# Put auxiliary/build files in the local temp folder instead of this synced folder.
# Only the final report.pdf (and its .synctex.gz) end up next to the sources.
use File::Spec;
my $tmp = $ENV{TEMP} || $ENV{TMPDIR} || '/tmp';
$aux_dir = File::Spec->catdir($tmp, 'latex-build', 'Thesis_report');
$out_dir = '.';
