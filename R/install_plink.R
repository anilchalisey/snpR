#' Function to download PLINK binary
#'
#' @return Path to PLINK binary
#' @export
#' @importFrom utils unzip download.file
#'
#' @examples
#' plink <- install_plink()

install_plink <- function() {
  
  os <- .Platform$OS.type
  if (os == "windows") {
    utils::download.file("https://www.cog-genomics.org/static/bin/plink180913/plink_win64.zip",
                         "plink.zip")
  } else {
    if (os == "darwin") {
      utils::download.file("https://www.cog-genomics.org/static/bin/plink180913/plink_mac.zip", 
                           "plink.zip")
    } else {
      utils::download.file("https://www.cog-genomics.org/static/bin/plink180913/plink_linux_x86_64.zip", 
                           "plink.zip")
    }
  }
  
  utils::unzip("plink.zip", exdir = "plink")
  return(file.path("plink", "plink.exe"))
}