# Keep generated artifacts out of the source tree so only manuscript inputs
# remain in `latex/` during normal editing.
$out_dir = 'build';
$aux_dir = 'build';

$pdf_mode = 1;
$bibtex_use = 2;
$max_repeat = 5;

$pdflatex = 'pdflatex -interaction=nonstopmode -file-line-error -synctex=1 %O %S';
$biber = 'biber --input-directory=build --output-directory=build %O %B';
