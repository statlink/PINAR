.onAttach <- function(lib, pkg) {
	version <- read.dcf(file.path(lib, pkg, "DESCRIPTION"),"Version")
	packageStartupMessage(paste("\nPNAR: ",version))
	msg <-  paste(
			r"( _ _ _ _    _     _    _ _ _ _    _ _ _ _ )",
			r"(|  _ _  |  | \   | |  |  ___  |  |  ___  |)",
			r"(| |_ _| |  |  \  | |  | |___| |  | |___| |)",
			r"(|  _ _ _|  |   \_| |  |  _ _  |  |    _ _|)",
			r"(| |        |  _    |  | |   | |  | |\ \   )",
			r"(| |        | | \   |  | |   | |  | | \ \_ )",
			r"(|_|        |_|  \__|  |_|   |_|  |_|  \__|)"
	,sep="\n")
	packageStartupMessage(msg)
}
