##### 20.11.2025 #####
install.packages("renv")
renv::init()       # sets up renv in your project
renv::snapshot()   # records the current package versions into renv.lock

#sink capture output into a text file
sink("notes/sessionInfo.txt")
print(sessionInfo())
sink()

#sink(file, append = TRUE) 
#if you want to add output to an existing file instead of overwriting it.
renv::deactivate() 
#I had problems with renv being on - could not install packages