# Shared by the vignettes: show cli's coloured console output (green/red decisions
# and intervals) in the rendered HTML. Only active when fansi is installed; otherwise
# the output is plain, so building the vignettes never depends on fansi.
# pkgdown converts the colour codes itself (and would undo this hook), so the hook is
# skipped on the website; the site's colours are set in pkgdown/extra.css.
if (!identical(Sys.getenv("IN_PKGDOWN"), "true") && requireNamespace("fansi", quietly = TRUE)) {
  options(cli.num_colors = 256L)
  # knitr's output hook (where the colours are converted) is bypassed when source and
  # output are collapsed into one block, so output goes into its own block
  knitr::opts_chunk$set(collapse = FALSE)
  # The eight basic terminal colours and their bright variants (solarized palette):
  # the same green and red as the asciicast "readme" theme used for the README.
  cols = c("#073642", "#DC322F", "#859900", "#B58900", "#268BD2", "#D33682", "#2AA198", "#EEE8D5",
           "#586E75", "#CB4B16", "#859900", "#B58900", "#268BD2", "#6C71C4", "#2AA198", "#FDF6E3")
  css = c("PRE.fansi SPAN {padding-top: .25em; padding-bottom: .25em};",
          sprintf(".fansi-color-%03d {color: %s;}", 0:15, cols))
  # cli prints the decision lines as messages, so hook messages and warnings as well
  fansi::set_knit_hooks(
    knitr::knit_hooks,
    which = c("output", "message", "warning"),
    proc.fun = function(x, class) {
      fansi::html_code_block(fansi::to_html(fansi::html_esc(x), classes = TRUE), class = class)
    },
    style = css
  )
}
