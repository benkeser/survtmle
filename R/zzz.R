.onAttach <- function(...) {
  packageStartupMessage("survtmle: Targeted Learning for Survival Analysis")
  packageStartupMessage("Version: ",
                        utils::packageDescription("survtmle")$Version)
}
